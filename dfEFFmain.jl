using LinearAlgebra
using Plots

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

T=24;
costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);


dfController=DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N,false,[1,2,3,4],[1,2,3,4]);
v,s,distSeq=runMCSimulation(dfController,demand,t,T,x0,costParam);

dfStdController=DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,false,[1,2,3,4],[1,2,3,4]);
v2,s2=run_WStdDfController(dfStdController,demand,t,T,x0,costParam,distSeq);
