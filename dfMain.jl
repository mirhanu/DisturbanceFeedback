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

function save_solution_time_plot(hors, times, title, filename; order=3)
    # Compute the polynomial fit coefficients
    hors_norm = hors / maximum(hors);
    coeffs = (hors_norm.^order) \ times;

    # Generate the plot
    p = plot(hors, 
        [times hors_norm.^order * coeffs], 
        xaxis=:log, yaxis=:log, 
        label=["Measured Time" "$(order)th Order Fit"], 
        xlabel="Horizon", ylabel="Per Iteration Time (s)", 
        title=title,
        legend=:topleft
    )
    
    # Save the plot
    savefig(p, filename)
end

# Efficient Controller Version
function compute_solution_times(controller::DFController, hors, demand; t=0, T=1, repeat=1)
    solveTimes = zeros(length(hors));
    perIterationTimes = zeros(length(hors));
    for i in 1:length(hors)
        for j in 1:repeat
            N = hors[i];
            set_N(controller,N);  # Set horizon

            # Run Efficient Controller simulation
            costParam.t = t;
            v, s, distSeq = runMCSimulation(controller, demand, t, T, x0, costParam);

            # Accumulate solve time
            solveTimes[i] += solve_time(controller.model);
            
            #Accumulate per iteration time (Note that barrier_iterations only works for Ipopt)
            sol_summary=solution_summary(controller.model)
            perIterationTimes[i] += solve_time(controller.model)/sol_summary.barrier_iterations;
        end

        # Average over repetitions
        solveTimes[i] /= repeat;
        perIterationTimes[i] /= repeat;
    end

    return solveTimes, perIterationTimes
end

# Standard Controller Version
function compute_solution_times(controller::DFStdController, hors, demand; t=0, T=1, repeat=1)
    solveTimes = zeros(length(hors));
    perIterationTimes = zeros(length(hors));

    for i in 1:length(hors)
        for j in 1:repeat
            N = hors[i];
            set_N(controller,N);  # Set horizon

            # Run Standard Controller simulation
            costParam.t = t;
            v2, s2 = runMCSimulation(controller, demand, t, T, x0, costParam, zeros(controller.n,1));

            # Accumulate solve time
            solveTimes[i] += solve_time(controller.model);

            #Accumulate per iteration time (Note that barrier_iterations only works for Ipopt)
            sol_summary=solution_summary(controller.model)
            perIterationTimes[i] += solve_time(controller.model)/sol_summary.barrier_iterations;
        end

        # Average over repetitions
        solveTimes[i] /= repeat;
    end

    return solveTimes, perIterationTimes
end


# ================================================================
# Run Monte Carlo simulations for both MPC formulations:
# - Efficient Formulation (DFController)
# - Standard Formulation (DFStdController)
#
# The results of both simulations should be identical.
# ================================================================

# Define the simulation horizon
# T = 2  

# # Set up cost parameters for the water distribution network
# costParam = WDNCostParams(sysC, sysD, elecPrice, reservoirPressures, t)

# # Run the Disturbance Feedback MPC with the Efficient Formulation
# dfController = DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N, false, [1,2,3,4], [1,2,3,4])
# inpsEff, statesEff, distSeq = runMCSimulation(dfController, demand, t, T, x0, costParam)

# # Run the Disturbance Feedback MPC with the Standard Formulation
# costParam.t = t  # Ensure the start time is the same as the Efficient one
# dfStdController = DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N, false, [1,2,3,4], [1,2,3,4])
# inpsSTD, statesSTD = runMCSimulation(dfStdController, demand, t, T, x0, costParam, distSeq)


# ================================================================
# Compute and plot the solution time for both MPC formulations 
# as the prediction horizon varies. This helps evaluate the 
# computational efficiency of the Efficient and Standard approaches.
# ================================================================

# # Define the range of prediction horizons to test
# hors=[4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60];

# # Set up cost parameters for the water distribution network
# costParam = WDNCostParams(sysC, sysD, elecPrice, reservoirPressures, t)

# # Small error bound to ensure feasibility for high prediction horizons
# E = [1e-7 0; 0 1e-7]  

# # Initialize the Efficient and Standard MPC controllers
# dfController = DFController(A, B1, B2[:,:], E, C, D, b, Y, z, N, false, [1,2,3,4], [1,2,3,4])
# dfStdController = DFStdController(A, B1, B2[:,:], E, C, D, b[:,:], Y, z[:,:], N, false, [1,2,3,4], [1,2,3,4])

# # Compute solution times for both controllers across different horizons
# effTim, effPerItTim = compute_solution_times(dfController, hors, demand; repeat=1)    # Efficient Controller
# stdTim, stdPerItTim = compute_solution_times(dfStdController, hors, demand; repeat=1) # Standard Controller

# # Generate and save the solution time plot for the Efficient Implementation (Cubic Fit)
# save_solution_time_plot(hors, effPerItTim, "Efficient Implementation", "efficient_implementation.png"; order=3)

# # Generate and save the solution time plot for the Standard Implementation (Sixth Order Fit)
# save_solution_time_plot(hors, stdPerItTim, "Standard Implementation", "standard_implementation.png"; order=6)