using JuMP
using Ipopt
using(LinearAlgebra)
using SparseArrays
using Plots
using OSQP
using Distributions
using Base.Threads
using JLD2


include("dFEffStruct.jl")
include("utils.jl")
include("parameters.jl")
include("stdMPCStruct.jl")
include("dfStdStruct.jl")
include("periodicMPC.jl")


 function runMCSimulation(controller::DFController,demand,t,T,x0,costParam,extreme=0)
    states=zeros(controller.n,T+1);
    inps=zeros(controller.m,T);
    states[:,1]=x0;
    disturbance=zeros(controller.n,1);
    distSeq=zeros(controller.n,T);
    v=0;
    s=0;
    for i=1:T
        costParam.t=t;
        demandinHorizon=get_ExInpinHorizon(controller.N,t,demand);
        v,s=calc_Inp(controller,states[:,i],demandinHorizon,costParam,disturbance);
        sample=rand(Uniform(-1,1), controller.l,1);
        if extreme==1
            sample=[-rand(Uniform(0.5,1));rand(Uniform(0.5,1))];
        elseif extreme==2
            sample=[-1;1];
        end
        disturbance=controller.E*sample;
        # sample=[-rand(Uniform(0.5,1));rand(Uniform(0.5,1))];
        # sample=[-1;1];
        # disturbance=controller.E*(sample);
        distSeq[:,i]=disturbance;
        states[:,i+1]=controller.A*states[:,i]+controller.B1*v[1:controller.m]+controller.B2*demandinHorizon[:,1]+disturbance;
        inps[:,i]=v[1:controller.m];
        t=t+1;
    end
    return inps,states,distSeq;
end


 function runMCSimulation(controller::mpcController,demand,t,T,x0,costParam,distSeq)
    states=zeros(controller.n,T+1);
    inps=zeros(controller.m,T);
    states[:,1]=x0;
    disturbance=zeros(controller.n,1);
    v=0
    s=0
    for i=1:T
        costParam.t=t;
        demandinHorizon=get_ExInpinHorizon(controller.N,t,demand);
        v,s=calc_Inp(controller,states[:,i],demandinHorizon,costParam,true,[1, 2, 3, 4],[1,2,3,4]);
        disturbance=distSeq[:,i];
        states[:,i+1]=controller.A*states[:,i]+controller.B1*v[1:controller.m]+controller.B2*demandinHorizon[:,1]+disturbance;
        inps[:,i]=v[1:controller.m];
        t=t+1;
    end
    return inps,states;
end

 function runMCSimulation(controller::DFStdController,demand,t,T,x0,costParam,distSeq)
    states=zeros(controller.n,T+1);
    inps=zeros(controller.m,T);
    states[:,1]=x0;
    disturbance=zeros(controller.n,1);
    v=0
    s=0
    for i=1:T
        costParam.t=t;
        demandinHorizon=get_ExInpinHorizon(controller.N,t,demand);
        v,s=calc_Inp(controller,states[:,i],demandinHorizon,costParam,disturbance);
        disturbance=distSeq[:,i];
        states[:,i+1]=controller.A*states[:,i]+controller.B1*v[1:controller.m]+controller.B2*demandinHorizon[:,1]+disturbance;
        inps[:,i]=v[1:controller.m];
        t=t+1;
    end
    return inps,states;
end

 function runMCSimulation(controller::periodicMPC,demand,t,T,x0,costParam,distSeq)
    states=zeros(controller.n,T+1);
    inps=zeros(controller.m,T);
    states[:,1]=x0;
    disturbance=zeros(controller.n,1);
    v=0
    s=0
    for i=1:T
        costParam.t=t;
        N=24-mod(t,24);
        set_Prediction_Horizon(controller,N);
        demandinHorizon=get_ExInpinHorizon(controller.N,t,demand);
        v,s=calc_Inp(controller,states[:,i],demandinHorizon,costParam,true,[1, 2, 3, 4],[1,2,3,4]);
        disturbance=distSeq[:,i];
        states[:,i+1]=controller.A*states[:,i]+controller.B1*v[1:controller.m]+controller.B2*demandinHorizon[:,1]+disturbance;
        inps[:,i]=v[1:controller.m];
        t=t+1;
    end
    return inps,states;
end

extreme=1
T=24*10;
# T=1;
costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);

days=10
vDist=zeros(days,2,T)
sDist=zeros(days,2,T+1)
costDist=zeros(days,1)
violDist=zeros(days,1)

vNom1=zeros(days,2,T)
sNom1=zeros(days,2,T+1)
costNom1=zeros(days,1)
violNom1=zeros(days,1)

vNom2=zeros(days,2,T)
sNom2=zeros(days,2,T+1)
costNom2=zeros(days,1)
violNom2=zeros(days,1)

vNom3=zeros(days,2,T)
sNom3=zeros(days,2,T+1)
costNom3=zeros(days,1)
violNom3=zeros(days,1)

vNom4=zeros(days,2,T)
sNom4=zeros(days,2,T+1)
costNom4=zeros(days,1)
violNom4=zeros(days,1)

vPer=zeros(days,2,T)
sPer=zeros(days,2,T+1)
costPer=zeros(days,1)
violPer=zeros(days,1)


Threads.@threads for i = 1:days
    x0=[rand(Uniform(x1low+0.4,x1high-0.4)),rand(Uniform(x2low+0.4,x2high-0.4))];
    costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
    dfController=DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N,true,[1,2,3,4],[1,2,3,4]);
    t1 = time();
    v,s,distSeq=runMCSimulation(dfController,demand,t,T,x0,costParam,extreme);
    t2 = time();
    vDist[i,:,:]=v
    sDist[i,:,:]=s

    costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
    normalController1=mpcController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], 48);
    t3 = time() ;
    v2,s2=runMCSimulation(normalController1,demand,t,T,x0,costParam,distSeq);
    t4 = time() ;
    vNom1[i,:,:]=v2
    sNom1[i,:,:]=s2

    costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
    normalController2=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+[-E[1,1];-E[1,1];-E[2,2];-E[2,2];0;0;0;0], Y, z[:,:], 48);
    t3 = time() ;
    v2,s2=runMCSimulation(normalController2,demand,t,T,x0,costParam,distSeq);
    t4 = time() ;
    vNom2[i,:,:]=v2
    sNom2[i,:,:]=s2

    costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
    normalController3=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+1.5*[-E[1,1];-E[1,1];-E[2,2];-E[2,2];0;0;0;0], Y, z[:,:], 48);
    t3 = time() ;
    v2,s2=runMCSimulation(normalController3,demand,t,T,x0,costParam,distSeq);
    t4 = time() ;
    vNom3[i,:,:]=v2
    sNom3[i,:,:]=s2

    costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
    normalController4=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+2*[-E[1,1];-E[1,1];-E[2,2];-E[2,2];0;0;0;0], Y, z[:,:], 48);
    t3 = time() ;
    v2,s2=runMCSimulation(normalController4,demand,t,T,x0,costParam,distSeq);
    t4 = time() ;
    vNom4[i,:,:]=v2
    sNom4[i,:,:]=s2

    costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
    periodicController=periodicMPC(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,xT);
    t3 = time() ;
    v3,s3=runMCSimulation(periodicController,demand,t,T,x0,costParam,distSeq);
    t4 = time() ;
    vPer[i,:,:]=v3
    sPer[i,:,:]=s3
    
end
costParam.t=t

for i=1:days
    costDist[i]=calc_cost(sDist[i,:,:],vDist[i,:,:],costParam);
    costNom1[i]=calc_cost(sNom1[i,:,:],vNom1[i,:,:],costParam);
    costNom2[i]=calc_cost(sNom2[i,:,:],vNom2[i,:,:],costParam)
    costNom3[i]=calc_cost(sNom3[i,:,:],vNom3[i,:,:],costParam)
    costNom4[i]=calc_cost(sNom4[i,:,:],vNom4[i,:,:],costParam)
    costPer[i]=calc_cost(sPer[i,:,:],vPer[i,:,:],costParam);
    violDist[i]=sum(sum((C*sDist[i,:,1:T]+D*vDist[i,:,:]-kron(ones(1,T),b)).>=0.0001,dims=1) .!=0);
    violNom1[i]=sum(sum((C*sNom1[i,:,1:T]+D*vNom1[i,:,:]-kron(ones(1,T),b)).>=0.0001,dims=1) .!=0);
    violNom2[i]=sum(sum((C*sNom2[i,:,1:T]+D*vNom2[i,:,:]-kron(ones(1,T),b)).>=0.0001,dims=1) .!=0);
    violNom3[i]=sum(sum((C*sNom3[i,:,1:T]+D*vNom3[i,:,:]-kron(ones(1,T),b)).>=0.0001,dims=1) .!=0);
    violNom4[i]=sum(sum((C*sNom4[i,:,1:T]+D*vNom4[i,:,:]-kron(ones(1,T),b)).>=0.0001,dims=1) .!=0);
    violPer[i]=sum(sum((C*sPer[i,:,1:T]+D*vPer[i,:,:]-kron(ones(1,T),b)).>=0.0001,dims=1) .!=0);
end
jldsave("sims2/ExtSim(0.5-1)wDemandSoft20percent.jld2"; T, days, vDist,sDist,vNom1,sNom1,vNom2,sNom2,vNom3,sNom3,vNom4,sNom4,vPer,sPer,violDist,violNom1,violNom2,violNom3,violNom4,violPer,costDist,costNom1,costNom2,costNom3,costNom4,costPer)



# using Distributed
# addprocs(6)
# for i=1:6
#     Threads.@spawn begin
#         costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
#         dfController=DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N,false,[1,2,3,4],[1,2,3,4]);
#         t1 = time();
#         v,s,distSeq=runMCSimulation(dfController,demand,t,T,x0,costParam);
#         t2 = time();
#         vDist[1,:,:]=v
#         sDist[1,:,:]=s
#     end

#     # dfStdController=DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,false,[1,2,3,4],[1,2,3,4]);
#     # v2,s2=run_WStdDfController(dfStdController,demand,t,T,x0,costParam,distSeq);

#     Threads.@spawn begin
#         costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
#         normalController1=mpcController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], 48);
#         t3 = time() ;
#         v2,s2=runMCSimulation(normalController1,demand,t,T,x0,costParam,distSeq);
#         t4 = time() ;
#         vNom1[1,:,:]=v2
#         sNom1[1,:,:]=s2
#     end

#     Threads.@spawn begin
#         costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
#         normalController2=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+[-0.054;-0.054;-0.083;-0.083;0;0;0;0], Y, z[:,:], 48);
#         t3 = time() ;
#         v2,s2=runMCSimulation(normalController2,demand,t,T,x0,costParam,distSeq);
#         t4 = time() ;
#         vNom2[1,:,:]=v2
#         sNom2[1,:,:]=s2
#     end

#     Threads.@spawn begin
#         costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
#         normalController3=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+1.5*[-0.054;-0.054;-0.083;-0.083;0;0;0;0], Y, z[:,:], 48);
#         t3 = time() ;
#         v2,s2=runMCSimulation(normalController3,demand,t,T,x0,costParam,distSeq);
#         t4 = time() ;
#         vNom3[1,:,:]=v2
#         sNom3[1,:,:]=s2
#     end

#     Threads.@spawn begin
#         costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
#         normalController4=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+2*[-0.054;-0.054;-0.083;-0.083;0;0;0;0], Y, z[:,:], 48);
#         t3 = time() ;
#         v2,s2=runMCSimulation(normalController4,demand,t,T,x0,costParam,distSeq);
#         t4 = time() ;
#         vNom4[1,:,:]=v2
#         sNom4[1,:,:]=s2
#     end

#     Threads.@spawn begin
#         costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
#         periodicController=periodicMPC(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,xT);
#         t3 = time() ;
#         v3,s3=runMCSimulation(periodicController,demand,t,T,x0,costParam,distSeq);
#         t4 = time() ;
#         vPer[1,:,:]=v3
#         sPer[1,:,:]=s3
#     end
# end
# costParam.t=t



# costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
# dfController=DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N,false,[1,2,3,4],[1,2,3,4]);
# t1 = time();
# v,s,distSeq=runMCSimulation(dfController,demand,t,T,x0,costParam);
# t2 = time() ;

# # dfStdController=DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,false,[1,2,3,4],[1,2,3,4]);
# # v2,s2=run_WStdDfController(dfStdController,demand,t,T,x0,costParam,distSeq);

# costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
# normalController=mpcController(A, B1, B2[:,:], E, C, D, b[:,:]+[-0.1;-0.1;-0.1;-0.1;0;0;0;0], Y, z[:,:], 48);
# t3 = time() ;
# v2,s2=runMCSimulation(normalController,demand,t,T,x0,costParam,distSeq);
# t4 = time() ;

# # costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
# # periodicController=periodicMPC(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,xT);
# # t3 = time() ;
# # v3,s3=runMCSimulation(periodicController,demand,t,T,x0,costParam,distSeq);
# # t4 = time() ;

# costParam.t=t


