#Computationally efficient implementation of the disturbance feedback method

using JuMP
using Ipopt
using LinearAlgebra
using SparseArrays
using OSQP

 mutable struct DFController
    #Struct properties
    A::Matrix{Float64}
    B1::Matrix{Float64}
    B2::Matrix{Float64}
    E::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    b::Vector{Float64}
    Y::Matrix{Float64}
    z::Vector{Float64}
    AA::SparseMatrixCSC{Float64, Int64}
    BB::SparseMatrixCSC{AffExpr, Int64}
    CC::SparseMatrixCSC{Float64, Int64}
    DD::SparseMatrixCSC{AffExpr, Int64}
    model::Model
    n::Int
    m::Int
    s::Int
    r::Int
    l::Int
    N::Int
    sc::Int #Multiplier of slack variablein the cost  function
    isSoft::Bool
    softInds::Vector{Float64}
    terminalSoftInds::Vector{Float64}
    mstruct::Int #Parameter defining the structure of the input policy
    lastM::Matrix{Float64} #M matrix in the last feasible solution
    lastv::Vector{Float64} #v vector in the last feasible solution
    lastDistSequence::Matrix{Float64} #Disturbance sequence since last feasible solution
    feasibilityTrack::Vector{Int} #Variable holding the feasibility of the past problems
    x::Matrix{AffExpr}
    soft::Vector{Any} #Slack variable
end
 
# function Definitions
 function DFController(A, B1, B2, E, C, D, b, Y, z, N,isSoft=false,softInds=[],terminalSoftInds=[],mstruct=0,sc=10000000)
  n=size(A,1);
  m=size(B1,2);
  s=size(C,1);
  r=size(Y,1);
  l=size(E,2);
  model = Model(Ipopt.Optimizer);
  # set_optimizer_attribute(model, "max_iter", 1000)
  obj = DFController(A, B1, B2, E, C, D, b, Y, z, [;;], [;;], [;;], [;;], model,
  n ,m, s, r, l, N,sc,isSoft,softInds,terminalSoftInds,mstruct,[;;],[],zeros(N*n,1),[],zeros(AffExpr,(1,n)),[]);
  set_Problem(obj);
  return obj
end

#Functions generating system matrices like A0,B0, etc.
 function gen_A0(dfC::DFController)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  A0=spzeros(N*n,N*(n+m));
  In=sparse(I,n,n);
  A0[1:n,1:n+m]=[dfC.B1 -In];
  for i=1:N-1
      A0[i*n+1:(i+1)*n,m+(m+n)*(i-1)+
      1:m+(m+n)*(i-1)+2*n+m]=[dfC.A dfC.B1 -In];
  end
  return A0;
end

 function gen_Ap(dfC::DFController,p::Int)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  l=dfC.l;
  s=dfC.s;
  r=dfC.r;
  k=Int.(floor((p-1)/l));
  In=sparse(I,n,n);
  if k<N-1
      Ap=zeros((N-k-1)*n,(N-1-k)*(m+s+n)+r);
      Ap[1:n,1:m+s+n]=[dfC.B1 zeros(n,s) -In];
      for i=1:dfC.N-k-2
          Ap[i*n+1:(i+1)*n,m+s+(n+m+s)*(i-1)+1:m+s+(n+m+s)*(i-1)+2*n+m+s]=
          [dfC.A dfC.B1 zeros(n,s) -In];
      end
  else
      Ap=[ In zeros(n,r) ];
  end
  return Ap;
end

 function gen_b0(dfC::DFController,x0,ExInp)
  N=dfC.N;
  n=dfC.n;
  A=dfC.A;
  B2=dfC.B2;
  b0=zeros(AffExpr,N*n,1);
  b0[1:n,1]=-A*x0-B2*ExInp[:,1];
  for i=1:N-1
      b0[i*n+1:(i+1)*n,1]=-B2*ExInp[:,i+1];
  end
  return b0;
end

 function gen_bp(dfC::DFController,p)
  N=dfC.N;
  n=dfC.n;
  A=dfC.A;
  E=dfC.E;
  l=dfC.l;
  k=Int.(floor((p-1)/l));
  j=mod(p-1,l)+1;
  if k<N-1
      bp=spzeros((N-k-1)*n,1);
      bp[1:n,1]=-A*E[:,j];
  else
      bp=E[:,j];
  end
  return bp;
end

 function gen_C0(dfC::DFController)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  s=dfC.s;
  r=dfC.r;
  C=dfC.C;
  D=dfC.D;
  Y=dfC.Y;
  C0=zeros(N*s+r,N*(n+m));
  In=sparse(I,n,n);
  C0[1:s,1:m]=D;
  for i=1:N-1
      C0[i*s+1:(i+1)*s,m+(m+n)*(i-1)+1:m+(m+n)*(i)]=[C D];
  end
  C0[N*s+1:N*s+r,N*(n+m)-n+1:N*(n+m)]=Y;
  return C0;
end

 function gen_Cp(dfC::DFController,p)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  s=dfC.s;
  r=dfC.r;
  l=dfC.l;
  C=dfC.C;
  D=dfC.D;
  Y=dfC.Y;
  k=Int.(floor((p-1)/l));
  Cbar=[C;-C];
  Dbar=[D;-D];
  Ybar=[Y;-Y];
  H=[-sparse(I,s,s);-sparse(I,s,s)];
  Hf=[-sparse(I,r,r);-sparse(I,r,r)];
  if k<N-1
      Cp=spzeros((N-k-1)*2*s+2*r,(N-1-k)*(m+s+n)+r); 
      Cp[1:2*s,1:m+s]=[Dbar H];
      for i=1:N-k-2
          Cp[i*2*s+1:(i+1)*2*s,m+s+(i-1)*(m+s+n)+1:m+s+(i)*(m+s+n)]=[Cbar Dbar H];
      end
      Cp[(N-k-1)*2*s+1:(N-k-1)*2*s+2*r,(N-1-k)*(m+s+n)-n+1:(N-1-k)*(m+s+n)+r]=[Ybar Hf];
  else
      Cp=[Ybar Hf];
  end
  return Cp;
end

 function gen_d0(dfC::DFController,x0)
  N=dfC.N;
  s=dfC.s;
  r=dfC.r;
  C=dfC.C;
  z=dfC.z;
  b=dfC.b;
  d0=spzeros(AffExpr,N*s+r,1);
  d0[1:s,1]=b-C*x0;
  for i=1:N-1
      d0[i*s+1:(i+1)*s,1]=b;
  end
  d0[N*s+1:N*s+r,1]=z;  
  return d0;
end

 function gen_dp(dfC::DFController,p)
  N=dfC.N;
  s=dfC.s;
  r=dfC.r;
  l=dfC.l;
  j=mod(p-1,dfC.l)+1;
  Cbar=[C;-C];
  k=Int.(floor((p-1)/l));
  dp=spzeros((N-k-1)*2*s+2*r,1);
  if k<N-1
    dp[1:2*dfC.s,1]=-Cbar*dfC.E[:,j];
  end
  return dp;
end

 function gen_Jp(dfC::DFController,p)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  s=dfC.s;
  r=dfC.r;
  l=dfC.l;
  k=Int.(floor((p-1)/l));
  if k<N-1
      Jp=spzeros(N*s+r,(N-1-k)*(m+s+n)+r);
      for i=0:N-1
          if i<=k
              continue;
          end
          ind=m+(i-(k+1))*(m+s+n);
          Jp[i*s+1:(i+1)*s,ind+1:ind+s]=sparse(I,s,s);
      end
      Jp[N*s+1:N*s+r,(N-1-k)*(m+s+n)+1:(N-1-k)*(m+s+n)+r]=sparse(I,r,r);
  else
      Jp=spzeros(N*s+r,r+n);
      Jp=[spzeros(N*s,r+n);spzeros(r,n) sparse(I,r,r)];
  end
  return Jp;
end

 function gen_bigA(dfC::DFController,xLength,bLength)
  Ai=spzeros(bLength,xLength);
  xInd=1;
  yInd=1;
  A0=gen_A0(dfC);
  xSize,ySize=size(A0);
  Ai[xInd:xInd+xSize-1,yInd:yInd+ySize-1]=A0;
  xInd=xInd+xSize;
  yInd=yInd+ySize;
  for p=1:dfC.N*dfC.l
    Ap=gen_Ap(dfC,p);
    xSize,ySize=size(Ap);
    Ai[xInd:xInd+xSize-1,yInd:yInd+ySize-1]=Ap;
    xInd=xInd+xSize;
    yInd=yInd+ySize;
  end
  dfC.AA=Ai;
end

 function gen_bigB(dfC::DFController,x0,ExInp,bLength)
  Bi=spzeros(AffExpr,bLength,1);
  xInd=1;
  B0=gen_b0(dfC,x0,ExInp);
  xSize=size(B0,1);
  Bi[xInd:xInd+xSize-1,1]=B0;
  xInd=xInd+xSize;
  for p=1:dfC.N*dfC.l
      Bp=gen_bp(dfC,p);
      xSize=size(Bp,1);
      Bi[xInd:xInd+xSize-1,1]=Bp;
      xInd=xInd+xSize;
  end
  dfC.BB=Bi;
end

 function gen_bigC(dfC::DFController,xLength,dLength)
  Ci=spzeros(dLength,xLength);
  xInd=1;
  yInd=1;
  C0=gen_C0(dfC);
  xSize,ySize=size(C0);
  c0xSize=xSize;
  Ci[xInd:xInd+xSize-1,yInd:yInd+ySize-1]=C0;
  xInd=xInd+xSize;
  yInd=yInd+ySize;
  for p=1:dfC.N*dfC.l
      Cp=gen_Cp(dfC,p);
      Jp=gen_Jp(dfC,p);
      xSize,ySize=size(Cp);
      Ci[xInd:xInd+xSize-1,yInd:yInd+ySize-1]=Cp;
      Ci[1:c0xSize,yInd:yInd+ySize-1]=Jp;
      xInd=xInd+xSize;
      yInd=yInd+ySize;
  end
  dfC.CC=Ci;
end

 function gen_bigD(dfC::DFController,x0,dLength)
  Di=spzeros(AffExpr,dLength,1);
  xInd=1;
  D0=gen_d0(dfC,x0);
  xSize=size(D0,1);
  Di[xInd:xInd+xSize-1,1]=D0;
  xInd=xInd+xSize;
  for p=1:dfC.N*dfC.l
      Dp=gen_dp(dfC,p);
      xSize=size(Dp,1);
      Di[xInd:xInd+xSize-1,1]=Dp;
      xInd=xInd+xSize;
  end
  dfC.DD=Di;
end

# function to get the δ\mathbf{c}
 function get_δ(dfC::DFController,x)
  bigJ=gen_Jp(dfC,1);
  for p=2:dfC.N*dfC.l
      Jp=gen_Jp(dfC,p);
      bigJ=[bigJ Jp];
  end
  return bigJ*x[dfC.N*(dfC.n+dfC.m)+1:length(x)];
end

# function to get \mathbf{U} matrix
 function get_u(dfC::DFController,x) 
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  s=dfC.s;
  r=dfC.r;
  l=dfC.l;
  bigU=zeros(eltype(x),m*N,l*N);
  start_ind=N*(n+m)+1;
  for p=1:N*l
      k=Int.(floor((p-1)/l));
      if k<N-1
          input_inds=(start_ind:start_ind+m-1)';
          input_inds=kron(ones(1,N-1-k),input_inds)+kron((0:m+n+s:(N-2-k)*(m+n+s))',ones(1,n));
          bigU[(k+1)*m+1:m*N,p]=x[Int.(input_inds)];
      end
      start_ind=start_ind+(N-1-k)*(m+s+n)+r;
  end
  return bigU;
end

 function get_Yp(dfC::DFController,p)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  s=dfC.s;
  r=dfC.r;
  l=dfC.l;
  C=dfC.C;
  D=dfC.D;
  Y=dfC.Y;
  k=Int.(floor((p-1)/l));
  if k<N-1
      Yp=zeros(N*s+r,(N-1-k)*(m+s+n)+n+r);
      for i=0:N-1
          if i<=k
              continue;
          end
          x_ind=(i-(k+1))*(m+s+n);
          inp_ind=n+(i-(k+1))*(m+s+n);
          Yp[i*s+1:(i+1)*s,x_ind+1:x_ind+n]=C;
          Yp[i*s+1:(i+1)*s,inp_ind+1:inp_ind+m]=D;
      end
      Yp[N*s+1:N*s+r,(N-1-k)*(m+s+n)+1:(N-1-k)*(m+s+n)+n]=Y;
  else
      Yp=[zeros(N*s,r+2*n);Y zeros(r,n) zeros(r,r)];
  end
  return Yp;
end

 function get_Y(dfC::DFController,x)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  l=dfC.l;
  j=mod(1-1,l)+1;
  bigY=get_Yp(dfC,1);
  xaug=[E[:,j];x[N*(n+m)+1:N*(n+m)+size(bigY,2)-n]];
  sz=size(bigY,2)-n;
  for p=2:N*l
      Yp=get_Yp(dfC,p);
      bigY=[bigY Yp];
      j=mod(p-1,l)+1;
      xaug=[xaug;E[:,j];x[N*(n+m)+sz+1:N*(n+m)+sz+size(Yp,2)-n]];
      sz=sz+size(Yp,2)-n;
  end
  return bigY*xaug;
end

# function to get the inputs and state predictions
 function get_StateInputs(dfC::DFController,x)
  N=dfC.N;
  m=dfC.m;
  n=dfC.n;
  state_inds=(m+1:m+n)';
  state_inds=kron(ones(1,N),state_inds)+kron((0:m+n:(N-1)*(m+n))',ones(1,n));
  input_inds=(1:m)';
  input_inds=kron(ones(1,N),input_inds)+kron((0:m+n:(N-1)*(m+n))',ones(1,n));
  states=vec(x[Int.(state_inds)]);
  v=vec(x[Int.(input_inds)]);
  return states,v;
end

 function calc_cost(dfC::DFController,states,Inps,x0,costParams)
  N=dfC.N;
  n=dfC.n;
  m=dfC.m;
  # Q=sparse(I,N*n,N*n);
  # R=sparse(I,N*m,N*m);
  # cost=Inps'*R*Inps+states'*Q*states;
  sysCb=kron(Matrix{Float64}(I, N, N),costParams.sysC);
  sysDb=kron(Matrix{Float64}(I, N, N),costParams.sysD);
  resPresB=kron(ones(N,1),costParams.reservoirPressures);
  heads=sysCb*[x0;states[1:(N-1)*n]]+sysDb*Inps;
  W=kron(Diagonal(vec(get_ExInpinHorizon(N,costParams.t,costParams.elecPrice))),Matrix{Float64}(I, n, n));
  cost=Inps'*W*(vec(heads)-vec(resPresB));
  if dfC.isSoft && length(dfC.softInds)>0
    cost=cost+dfC.sc*sum(dfC.soft);
  end
  return cost;
end

# function to put the slack variables in vector form.
 function construct_softVariable(dfC::DFController,liftedSoftInds)
  N=dfC.N;
  s=dfC.s;
  r=dfC.r;
  softLifted=zeros(AffExpr,N*s+r,1);
  softCount=1;
  for i=1:length(liftedSoftInds)
      softLifted[Int.(liftedSoftInds[i])]=dfC.soft[Int.(softCount)];
      softCount=softCount+1;
  end
  for i=1:length(dfC.terminalSoftInds)
    softLifted[N*s+Int.(dfC.terminalSoftInds[i])]=dfC.soft[Int.(softCount)];
    softCount=softCount+1;
  end
  return softLifted;
end
 function get_matrixSizes(dfC::DFController)
  N=dfC.N;
  l=dfC.l;
  n=dfC.n;
  m=dfC.m;
  s=dfC.s;
  r=dfC.r;
  xLength=0;
  x0Length=N*(n+m);
  xLength+=x0Length;
  bLength=0;
  b0Length=N*n;
  bLength+=b0Length;
  dLength=0;
  d0Length=N*s+r;
  dLength+=d0Length;
  for p=1:N*l
    k=Int.(floor((p-1)/l));
    dLength+=(N-k-1)*2*s+2*r;
    if k<N-1
      xLength+=(N-1-k)*(m+s+n)+r;
      bLength+=(N-k-1)*n;
    else
      xLength+=n+r;
      bLength+=size(dfC.E,1);
    end
  end
  return xLength,bLength,dLength
end
# function to set the optimization problem. Should be called anytime system parameters are changed.
 function set_Problem(dfC::DFController)
  #Initial condition x0 parameter
  @variable(dfC.model, Px0[i = 1:dfC.n] in Parameter(i));
  nExInp=size(dfC.B2,2);
  #External Input parameter
  @variable(dfC.model, PExInp[1:nExInp,1:dfC.N] in Parameter(0.0));
  #Generate matrices
  xLength,bLength,dLength=get_matrixSizes(dfC);
  gen_bigA(dfC,xLength,bLength);
  gen_bigB(dfC,Px0,PExInp,bLength);
  gen_bigC(dfC,xLength,dLength);
  gen_bigD(dfC,Px0,dLength);
  #Set constraints
  @variable(dfC.model,pre_x[1:Int.(size(dfC.AA,2)-dfC.l*dfC.m*(dfC.mstruct*(2*dfC.N-1-dfC.mstruct)/2))]);
  dfC.x=construct_x(dfC,pre_x);
  @constraint(dfC.model, c1, dfC.AA*dfC.x.==dfC.BB);
  if dfC.isSoft && length(dfC.softInds)>0
    liftedSoftInds=kron(ones(1,dfC.N),(dfC.softInds)')+kron((0*dfC.s:dfC.s:(dfC.N-1)*dfC.s)',ones(1,length(dfC.softInds)));
    @variable(dfC.model,soft[1:length(liftedSoftInds)+length(dfC.terminalSoftInds)]);
    dfC.soft=soft;
    softLifted=construct_softVariable(dfC,liftedSoftInds);
    Dsoft=dfC.DD[1:dfC.N*dfC.s+dfC.r,1]+softLifted;
    Dsoft=[Dsoft; convert(Vector{AffExpr},dfC.DD[dfC.N*dfC.s+dfC.r+1:length(dfC.DD)])];
    @constraint(dfC.model, c2, dfC.CC*dfC.x.<=Dsoft);
    @constraint(dfC.model,cSoft, dfC.soft.>=0)
  else
    @constraint(dfC.model, c2, dfC.CC*dfC.x.<=dfC.DD)
  end
end

# function to put x variable in vector form
 function construct_x(dfC::DFController,pre_x)
  N=dfC.N;
  l=dfC.l;
  m=dfC.m;
  n=dfC.n;
  s=dfC.s;
  r=dfC.r;
  mstruct=dfC.mstruct;
  x=zeros(AffExpr,(size(dfC.AA,2),1));
  x[1:N*(n+m)]=pre_x[1:N*(n+m)];
  pre_ind=N*(n+m)+1;
  x_ind=N*(n+m)+1;
  for p=1:N*l
    k=Int.(floor((p-1)/l));
    if k<N-1
      len_xp=(N-1-k)*(m+s+n)+r;
      for i=1:len_xp
        if mod(i,m+n+s)<=m && mod(i,m+n+s)>0 && len_xp-i>r && floor(i/(m+n+s))<mstruct
          x[x_ind]=0;
        else
          x[x_ind]=pre_x[pre_ind];
          pre_ind=pre_ind+1;
        end
        x_ind=x_ind+1;
      end
    else
      len_xp=n+r;
      for i=1:len_xp
        x[x_ind]=pre_x[pre_ind];
        pre_ind=pre_ind+1;
        x_ind=x_ind+1;
      end
    end
  end
  return x;
end

 function calc_Inp(dfC::DFController,x0,ExInp,costParams,disturbance)
  #Set BB and DD according to x0 and ExInp
  x=dfC.x;
  # dfC.BB[1:dfC.n*dfC.N,1]=gen_b0(dfC,x0,ExInp);
  # dfC.DD[1:dfC.N*dfC.s+dfC.r,1]=gen_d0(dfC,x0);

  states,v=get_StateInputs(dfC,x);
  cost=calc_cost(dfC,states,v,x0,costParams);
  # U=get_u(x);
  set_silent(dfC.model)

  set_parameter_value.(dfC.model[:Px0], x0);
  
  set_parameter_value.(dfC.model[:PExInp], ExInp);
  @objective(dfC.model, Min, cost);
  optimize!(dfC.model)
  if is_solved_and_feasible(dfC.model)
    push!(dfC.feasibilityTrack,1);
    U=get_u(dfC,JuMP.value.(x));
    J=kron(Matrix{Float64}(I,dfC.N,dfC.N),dfC.E);
    dfC.lastM=U*(J\ Matrix{Float64}(I,N*n,N*n));
    dfC.lastv= JuMP.value.(v);
    dfC.lastDistSequence=zeros(dfC.N*dfC.n,1);
    return JuMP.value.(v),JuMP.value.(states);
  else
    nSteps=length(dfC.feasibilityTrack);
    #If the problem is infeasible
    distInd=1;
    for i=1:nSteps
      if dfC.feasibilityTrack[nSteps-i+1]==1
        distInd=i;
        break;
      end
    end
    #Set the last disturbance
    dfC.lastDistSequence[(distInd-1)*dfC.n+1:distInd*dfC.n]=disturbance;
    push!(dfC.feasibilityTrack,0);
    #Return the input according to input policy found in the last feasible step
    Inp=dfC.lastM*dfC.lastDistSequence+dfC.lastv;
    currentInP=Inp[(distInd)*dfC.m+1:(distInd+1)*dfC.m];
    return currentInP, [dfC.A*x0+dfC.B1*currentInP+dfC.B2*ExInp[:,1]]
  end
end