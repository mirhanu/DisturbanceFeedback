using LinearAlgebra
using Plots
using Distributions

include("dFEffStruct.jl")
include("utils.jl")
include("parameters.jl")
include("dfStdStruct.jl")


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

#Simulation horizon
T=2;
costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);

#Run the Disturbance Feedback MPC with Efficient Formulation
dfController=DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N,false,[1,2,3,4],[1,2,3,4]);
inpsEff,statesEff,distSeq=runMCSimulation(dfController,demand,t,T,x0,costParam);

#Run the Disturbance Feedback MPC with Standard Formulation 
# (The results of both simulations should be the same.)
costParam.t=t;t;
dfStdController=DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,false,[1,2,3,4],[1,2,3,4]);
inpsSTD,statesSTD=runMCSimulation(dfStdController,demand,t,T,x0,costParam,distSeq);


# Uncomment the code below to plot the computation time for both 
# formulations as the prediction horizon varies
# hors=[4 8 12 16 20 24 28 32 36 40 44 48 52 56 60];
# costParam=WDNCostParams(sysC,sysD,elecPrice,reservoirPressures,t);
# dfController=DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N,false,[1,2,3,4],[1,2,3,4]);
# dfStdController=DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N,false,[1,2,3,4],[1,2,3,4]);
# effTim=zeros(1,length(hors));
# stdTim=zeros(1,length(hors));
# T=24
# repeat=1;
# for i=1:length(hors)
#     N=hors[i]
#     for j=1:repeat
#         costParam.t=t;
#         v,s,distSeq=runMCSimulation(dfController,demand,t,T,x0,costParam);
#         costParam.t=t;
#         v2,s2=runMCSimulation(dfStdController,demand,t,T,x0,costParam,distSeq);
#         effTim[i]+= solve_time(dfController.model);
#         stdTim[i]+= solve_time(dfStdController.model);
#     end
#     effTim[i]=effTim[i]/repeat;
#     stdTim[i]=stdTim[i]/repeat;
# end