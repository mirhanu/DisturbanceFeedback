#The standart implementation of the disturbance feedback method

using JuMP
using Ipopt
using LinearAlgebra
using SparseArrays
using OSQP

mutable struct DFStdController
    A::Matrix{Float64}
    B1::Matrix{Float64}
    B2::Matrix{Float64}
    E::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    b::Matrix{Float64}
    Y::Matrix{Float64}
    z::Matrix{Float64}
    Ab::SparseMatrixCSC{Float64, Int64}
    B1b::SparseMatrixCSC{Float64, Int64}
    B2b::SparseMatrixCSC{Float64, Int64}
    Cb::SparseMatrixCSC{Float64, Int64}
    Db::SparseMatrixCSC{Float64, Int64}
    F::SparseMatrixCSC{Float64, Int64}
    G::SparseMatrixCSC{Float64, Int64}
    T::SparseMatrixCSC{Float64, Int64}
    csmall::SparseMatrixCSC{Float64, Int64}
    J::SparseMatrixCSC{Float64, Int64}
    model::Model
    n::Int
    m::Int
    s::Int
    r::Int
    l::Int
    N::Int
    sc::Int #Multiplier of slack variable in the cost function
    isSoft::Bool
    softInds::Vector{Float64}
    terminalSoftInds::Vector{Float64}
    mstruct::Int #Parameter defining the structure of the input policy
    lastM::Matrix{Float64} #M matrix in the last feasible solution
    lastv::Vector{Float64} #v vector in the last feasible solution
    lastDistSequence::Matrix{Float64} #Disturbance sequence since last feasible solution
    feasibilityTrack::Vector{Int} #Variable holding the feasibility of the past problems
    states::Matrix{AffExpr}
    U::Matrix{AffExpr}
    λ::Matrix{AffExpr}
    v::Vector{AffExpr}
    soft::Vector{Any} #Slack variable
end 

#Function Definitions
function DFStdController(A, B1, B2, E, C, D, b, Y, z, N,isSoft=false,softInds=[],terminalSoftInds=[],mstruct=0,sC=10000000)
    n=size(A,1);
    m=size(B1,2);
    s=size(C,1);
    r=size(Y,1);
    l=size(E,2);
    model = Model(Ipopt.Optimizer);
    obj = DFStdController(A, B1, B2, E, C, D, b, Y, z, [;;],[;;],
     [;;], [;;], [;;], [;;], [;;], [;;], [;;], [;;], model,
     n ,m, s, r, l, N,sC,isSoft,softInds,terminalSoftInds,mstruct,[;;],[],spzeros(N*n,1),[],spzeros(AffExpr,(1,n)),spzeros(AffExpr,(1,n)),spzeros(AffExpr,(1,n)),[],[]);
    set_Problem(obj);
    return obj
end

#Set the prediction horizon
function set_N(dfStdC::DFStdController,N)
	dfStdC.N=N;
	set_Problem(dfStdC);
end

#Functions generating system matrices like A0,B0, etc.
function gen_bigA(dfStdC::DFStdController)
    n=dfStdC.n;
    N=dfStdC.N;
    A=dfStdC.A;
    Ab=spzeros((N+1)*n, n);
    ai=sparse(I,n,n);
    for i in 1:n:n*(N+1)-1
        Ab[i:i+n-1,:]=ai;
        ai=A*ai;
    end
    dfStdC.Ab=Ab;
end

function gen_bigE(dfStdC::DFStdController)
    n=dfStdC.n;
    N=dfStdC.N;
    Ab=dfStdC.Ab;
    Su=spzeros((N+1)*n,n*N);
    columns=Ab[1:n*N,:];
    for i=1:n:n*N-n+1
        Su[i+n:n*(N+1),i:i+n-1]=columns[1:n*N-i+1,:];
    end
    return Su;
end

function set_bigC(dfStdC::DFStdController)
    n=dfStdC.n;
    N=dfStdC.N;
    Y=dfStdC.Y;
    C=dfStdC.C;
    dfStdC.Cb=[kron(sparse(I,N,N),C) spzeros(N*size(C,1),n);spzeros(size(Y,1),n*N) Y];
end
function set_bigD(dfStdC::DFStdController)
    N=dfStdC.N;
    m=dfStdC.m;
    Y=dfStdC.Y;
    D=dfStdC.D;
    dfStdC.Db=[kron(sparse(I,N,N),D); spzeros(size(Y,1),m*N)];
end
function set_bigB1(dfStdC::DFStdController,Eb)
    N=dfStdC.N;
    B1=dfStdC.B1;
    dfStdC.B1b=Eb*kron(sparse(I,N,N),B1);
end
function set_bigB2(dfStdC::DFStdController,Eb)
    N=dfStdC.N;
    B2=dfStdC.B2;
    dfStdC.B2b=Eb*kron(sparse(I,N,N),B2);
end
function set_csmall(dfStdC::DFStdController)
    N=dfStdC.N;
    b=dfStdC.b;
    z=dfStdC.z;
    dfStdC.csmall=[kron(ones(N,1),b);z];
end
function set_F(dfStdC::DFStdController)
    Cb=dfStdC.Cb;
    B1b=dfStdC.B1b;
    Db=dfStdC.Db;
    dfStdC.F=Cb*B1b+Db;
end
function set_T(dfStdC::DFStdController)
    Cb=dfStdC.Cb;
    Ab=dfStdC.Ab;
    dfStdC.T=-Cb*Ab;
end
function set_G(dfStdC::DFStdController,Eb)
    Cb=dfStdC.Cb;
    dfStdC.G=Cb*Eb;
end
function set_J(dfStdC::DFStdController)
    E=dfStdC.E;
    N=dfStdC.N;
    dfStdC.J=kron(sparse(I,N,N),E);
end

function set_U(dfStdC::DFStdController,uij)
    n=dfStdC.n;
    m=dfStdC.m;
    l=dfStdC.l;
    N=dfStdC.N;
    mstruct=dfStdC.mstruct;
    U=zeros(AffExpr,m*N,l*N);
    uCount=1;
    for i=0:N-2
        startInd=mstruct+i+1;
        if startInd>=N
            break;
        end
        for j=1:l
            U[startInd*m+1:m*N,i*l+j]=uij[uCount:uCount+m*(N-startInd)-1];
            uCount=uCount+m*(N-startInd);
        end
    end
    return U;
end

#Function to set the optimization problem.
function set_Problem(dfStdC::DFStdController)
    empty!(dfStdC.model);
    #Generate matrices
    gen_bigA(dfStdC);
    Eb=gen_bigE(dfStdC);
    set_bigC(dfStdC);
    set_bigD(dfStdC);
    set_bigB1(dfStdC,Eb);
    set_bigB2(dfStdC,Eb);
    set_csmall(dfStdC);
    set_F(dfStdC);
    set_T(dfStdC);
    set_G(dfStdC,Eb);
    set_J(dfStdC);
    
    s=dfStdC.s;
    n=dfStdC.m;
    m=dfStdC.m;
    l=dfStdC.l;
    N=dfStdC.N;
    r=dfStdC.r;
    mstruct=dfStdC.mstruct;
    Ab=dfStdC.Ab;
    B1b=dfStdC.B1b;
    B2b=dfStdC.B2b;
    Cb=dfStdC.Cb;
    F=dfStdC.F;
    J=dfStdC.J;
    G=dfStdC.G;
    csmall=dfStdC.csmall;
    T=dfStdC.T;
    #Optimization variables
    @variable(dfStdC.model, uij[1:Int.(m*l*(N*(N-1)/2-mstruct*(2*N-1-mstruct)/2))]);
    U=set_U(dfStdC,uij);
    # @variable(dfStdC.model, U[1:m*N,1:l*N]);
    @variable(dfStdC.model, λ[1:s*N+r,1:l*N]);
    @variable(dfStdC.model, v[1:m*N]);
    dfStdC.U=U;
    dfStdC.λ=λ;
    dfStdC.v=v;
    #Initial condition x0 parameter
    @variable(dfStdC.model, Px0[i = 1:n] in Parameter(i));
    nExInp=size(dfStdC.B2,2);
    #External Input parameter
    @variable(dfStdC.model, PExInp[1:nExInp,1:N] in Parameter(0.0));
    dfStdC.states=Ab*Px0+B1b*v+B2b*PExInp';
    #Set constraints
    if dfStdC.isSoft && length(dfStdC.softInds)>0
        liftedSoftInds=kron(ones(1,N),(dfStdC.softInds)')+kron((0:s:(N-1)*s)',ones(1,length(dfStdC.softInds)));
        @variable(dfStdC.model,soft[1:length(liftedSoftInds)+length(dfStdC.terminalSoftInds)]);
        dfStdC.soft=soft;
        softLifted=construct_softVariable(dfStdC,liftedSoftInds)
        @constraint(dfStdC.model, c1, F*v+λ*ones(l*N,1) .<=  csmall+T*x0-Cb*B2b*PExInp'+softLifted);
        @constraint(dfStdC.model,cSoft, soft.>=0);
    else
        @constraint(dfStdC.model, c1, F*v+λ*ones(l*N,1) .<=  csmall+T*Px0-Cb*B2b*PExInp')
    end
    @constraint(dfStdC.model, c2, F*U+G*J .>= -λ)
    @constraint(dfStdC.model, c3, F*U+G*J .<=  λ);
end

#Function to put the slack variables in vector form.
function construct_softVariable(dfStdC::DFStdController,liftedSoftInds)
    N=dfStdC.N;
    s=dfStdC.s;
    r=dfStdC.r;
    softLifted=spzeros(AffExpr,N*s+r,1);
    softCount=1;
    for i=1:length(liftedSoftInds)
        softLifted[Int.(liftedSoftInds[i])]=dfStdC.soft[Int.(softCount)];
        softCount=softCount+1;
    end
    for i=1:length(dfStdC.terminalSoftInds)
        softLifted[N*s+Int.(dfStdC.terminalSoftInds[i])]=dfStdC.soft[Int.(softCount)];
        softCount=softCount+1;
    end
    return softLifted;
end

#Function to set constraints on U matrix
function uConst(dfStdC::DFStdController,U)
    l=dfStdC.l;
    m=dfStdC.m;
    row=size(U,1);
    col=size(U,2);
    for i=1:row
        for j=1:col
            inpNo=floor((i-1)/m)+1;
            distNo=floor((j-1)/l)+1;
            if inpNo<=distNo+dfStdC.mstruct
                @constraint(dfStdC.model,  U[i,j]==0);
            end
        end
    end
end
function calc_cost(dfStdC::DFStdController,states,Inps,costParams)
    N=dfStdC.N;
    n=dfStdC.n;
    sysCb=kron(Matrix{Float64}(I, N, N),costParams.sysC);
    sysDb=kron(Matrix{Float64}(I, N, N),costParams.sysD);
    resPresB=kron(ones(N,1),costParams.reservoirPressures);
    heads=sysCb*states[1:N*n,1]+sysDb*Inps;
    W=kron(Diagonal(vec(get_ExInpinHorizon(N,costParams.t,costParams.elecPrice))),Matrix{Float64}(I, n, n));
    cost=Inps'*W*(vec(heads)-vec(resPresB));
    if dfStdC.isSoft && length(dfStdC.softInds)>0
        cost=cost+dfStdC.sc*sum(dfStdC.soft);
    end
    return cost;
end
function calc_Inp(dfStdC::DFStdController,x0,ExInp,costParams,disturbance)
    N=dfStdC.N;
    m=dfStdC.m;
    n=dfStdC.n;
    J=dfStdC.J;

    cost=calc_cost(dfStdC,dfStdC.states,dfStdC.v,costParams);
    set_parameter_value.(dfStdC.model[:Px0], x0);
    set_parameter_value.(dfStdC.model[:PExInp], ExInp);

    @objective(dfStdC.model, Min, cost)
    optimize!(dfStdC.model)
    if is_solved_and_feasible(dfStdC.model)
        push!(dfStdC.feasibilityTrack,1);
        U=JuMP.value.(dfStdC.U);
        dfStdC.lastM=U*(Matrix(J)\ Matrix{Float64}(I,N*n,N*n));
        dfStdC.lastv= JuMP.value.(dfStdC.v);
        dfStdC.lastDistSequence=spzeros(N*n,1);
        return JuMP.value.(dfStdC.v),JuMP.value.(dfStdC.states);
    else
        nSteps=length(dfStdC.feasibilityTrack);
        #If the problem is infeasible
        distInd=1;
        for i=1:nSteps
            if dfStdC.feasibilityTrack[nSteps-i+1]==1
            distInd=i;
            break;
            end
        end
        #Set the last disturbance
        dfStdC.lastDistSequence[(distInd-1)*n+1:distInd*n]=disturbance;
        push!(dfStdC.feasibilityTrack,0);
        #Return the input according to input policy found in the last feasible step
        Inp=dfStdC.lastM*dfStdC.lastDistSequence+dfStdC.lastv;
        currentInP=Inp[(distInd)*m+1:(distInd+1)*m];
        return currentInP, [dfStdC.A*x0+dfStdC.B1*currentInP+dfStdC.B2*ExInp[:,1]]
    end
end
  