module Simulate

using DifferentialEquations
using Dates
using DataFrames
using DataInterpolations

# Import siblings
using ...Shared.TimeUtil  # <--- Now uses Global Time
using ..Constants
using ..MosquitoModelDynamics

export run_mosquito_simulation

"""
    run_mosquito_simulation(trap_df, temp_interp)

Runs a forward simulation using Global Time (t=0 is GLOBAL_START_DATE).
"""
function run_mosquito_simulation(trap_df::DataFrame,
                                 temp_interp;
                                 solver = Rosenbrock23(), # Changed default to robust solver
                                 saveat_daily::Float64 = 1.0,
                                 reltol::Float64 = 1e-6,
                                 abstol::Float64 = 1e-6)

    println("\n=== DIAGNOSTIC: Global Time Alignment ===")

    # 1. Setup Time Axis (Using Global t)
    t_start_data = date_to_t(minimum(trap_df.date))
    t_end_data   = date_to_t(maximum(trap_df.date))

    println("1. Global Epoch:        $GLOBAL_START_DATE")
    println("2. Data Range (t):      $t_start_data to $t_end_data")
    
    # Check Temp alignment (Is it summer temp in summer?)
    temp_check = temp_interp(t_start_data)
    println("3. Temp at start (t=$t_start_data): $(round(temp_check, digits=2))°C")

    # 2. Set Simulation Span 
    # Buffer: Start 30 days early to let populations stabilize
    tspan = (t_start_data - 30.0, t_end_data + 10.0)

    # 3. Construct Parameters
    params = MosquitoModelDynamics.MosquitoModelParams(temp_interp = temp_interp)
    u0 = MosquitoModelDynamics.default_u0(params)
    
    prob = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0, tspan, params)

    println("--- Starting Solver ---")
    sol = solve(prob, solver; saveat=saveat_daily, reltol=reltol, abstol=abstol)

    if Symbol(sol.retcode) != :Success
        @warn "Simulation solver failed with status: $(sol.retcode)"
    end

    # 4. Generate Outputs

    # A. Smooth Trajectories
    daily_A = [u[MosquitoModelDynamics.IX_A] for u in sol.u]
    daily_M = [u[MosquitoModelDynamics.IX_M] for u in sol.u]

    # Map global t back to real Dates
    sim_dates = t_to_date.(sol.t)

    # B. Discrete Predictions
    # Get the exact global times for the trap data points
    data_times = date_to_t.(trap_df.date)

    mfai_pred = MosquitoModelDynamics.compute_mfai_theo(
        sol,
        data_times,
        MosquitoModelDynamics.IX_T,
        Constants.N_TRAPS
    )

    println("Simulation complete. Steps: $(length(sol.t))")

    return (
        times = sol.t,
        dates = sim_dates,
        A = daily_A,
        M = daily_M,
        mfai_times = data_times,
        mfai_pred = mfai_pred,
        sol_obj = sol,
        params = params
    )
end

end # module


# module Simulate

# using DifferentialEquations
# using Dates
# using DataFrames

# # Import siblings
# using ..Constants
# using ..MosquitoModelDynamics

# export run_mosquito_simulation

# """
#     run_mosquito_simulation(trap_df, temp_interp, fitted_params)

# Runs a forward simulation of the mosquito model using the optimized parameters.

# Arguments:
# - `trap_df`: The processed DataFrame containing the trap data dates.
# - `temp_interp`: The temperature interpolation object.
# - `fitted_params`: Vector [C₀, bₖ, ϵ] resulting from the fitting process.

# Returns:
# - `simulation_results`: A NamedTuple containing:
#     - `times`: Vector of days (0 to end).
#     - `dates`: Vector of Date objects corresponding to `times`.
#     - `A`: Vector of Aquatic population (daily).
#     - `M`: Vector of Adult Female population (daily).
#     - `mfai_pred`: Vector of predicted MFAI values matching `trap_df` rows.
# """
# function run_mosquito_simulation(trap_df::DataFrame, 
#                                  temp_interp, 
#                                  fitted_params::Vector{Float64};
#                                  # CHANGE: expose solver + tolerances + saveat as kwargs
#                                  solver = Tsit5(),
#                                  saveat_daily::Float64 = 1.0,
#                                  reltol::Float64 = 1e-5,
#                                  abstol::Float64 = 1e-5)
    
#     println("--- Running Forward Simulation ---")

#     # 1. Setup Time Axis
#     start_date = minimum(trap_df.date)
#     end_date   = maximum(trap_df.date)
#     total_days = Dates.value(end_date - start_date)
    
#     # We simulate a few days past the end to ensure plots don't cut off
#     tspan = (0.0, Float64(total_days + 5))

#     # 2. Reconstruct Parameters
#     # We use the fitted values [C0, bk, eps] + fixed constants
#     C0, bk, eps = fitted_params
    
#     params = MosquitoModelDynamics.MosquitoModelParams(
#         Nⱼ = Constants.N_TRAPS,
#         H₀ = Constants.N_HOUSEHOLDS,
#         j  = Constants.ALPHA,
#         k  = Constants.K,
#         C₀ = C0,
#         bₖ = bk,
#         ϵ  = eps,
#         temp_interp = temp_interp
#     )

#     # 3. Setup and Solve ODE
#     # Note: Initial conditions A(0) depend on C0, so we recalculate them
#     u0 = MosquitoModelDynamics.default_u0(params)
#     prob = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0, tspan, params)

#     # Solve with high resolution (daily) for smooth plotting
#     # saveat=1.0 gives us the state for every single day
#     # CHANGE: use kwargs saveat_daily/reltol/abstol/solver instead of hardcoding
#     sol = solve(prob, solver, saveat=saveat_daily, reltol=reltol, abstol=abstol)

#     if sol.retcode != :Success
#         @warn "Simulation solver failed with status: $(sol.retcode)"
#     end

#     # 4. Generate Outputs
    
#     # A. Smooth Trajectories (for Line Plots)
#     # sol.t is the time in days, sol[i] is the state vector at that time
#     daily_A = [u[MosquitoModelDynamics.IX_A] for u in sol.u]
#     daily_M = [u[MosquitoModelDynamics.IX_M] for u in sol.u]
    
#     # Map simulation days back to real Dates
#     # CHANGE: use round(Int, t) instead of floor(t) to avoid FP off-by-one issues
#     sim_dates = [start_date + Day(round(Int, t)) for t in sol.t]

#     # B. Discrete Predictions (for Scatter Plots vs Data)
#     # We need the model's predicted MFAI specifically at the data points
#     data_times = Float64[Dates.value(d - start_date) for d in trap_df.date]
    
#     # Use the helper function from ModelDynamics to ensure math is identical to fitting
#     mfai_pred = MosquitoModelDynamics.compute_mfai_theo(
#         sol, 
#         data_times, 
#         MosquitoModelDynamics.IX_T, 
#         Constants.N_TRAPS
#     )

#     println("Simulation complete. Generated $(length(daily_M)) daily points.")

#     return (
#         times = sol.t,
#         dates = sim_dates,
#         A = daily_A,
#         M = daily_M,
#         # CHANGE: return the data_times too so Report doesn't recompute it
#         mfai_times = data_times,
#         mfai_pred = mfai_pred,
#         sol_obj = sol # Return full solution just in case
#     )
# end

# end # module
