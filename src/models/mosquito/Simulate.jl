module Simulate

using DifferentialEquations
using Dates
using DataFrames

using ..Constants
using ..MosquitoModelDynamics

export run_mosquito_simulation

"""
    run_mosquito_simulation(trap_df, temp_interp; fitted_params=nothing, kwargs...)

Runs a forward simulation of the mosquito model.
automatically aligns t=0 with the earliest date in trap_df.

Arguments:
- `trap_df`: DataFrame containing :date column.
- `temp_interp`: Temperature interpolation function (t -> T).
- `fitted_params`: Optional [C₀, bₖ, ϵ]. If nothing, uses defaults.

Returns a NamedTuple with aligned times, dates, state trajectories, and trap predictions.
"""
function run_mosquito_simulation(trap_df::DataFrame,
                                 temp_interp;
                                 fitted_params::Union{Nothing, Vector{Float64}} = nothing,
                                 solver = Tsit5(),
                                 saveat_daily::Float64 = 1.0,
                                 reltol::Float64 = 1e-6,
                                 abstol::Float64 = 1e-6)

    println("\n=== Starting Simulation ===")

    # 1. Setup Time Axis (Align with Fitting.jl logic)
    #    t=0 corresponds to the first day in the data.
    start_date = minimum(trap_df.date)
    end_date   = maximum(trap_df.date)
    total_days = Float64(Dates.value(end_date - start_date))

    println("Simulation Range: $start_date to $end_date ($total_days days)")

    # Simulate 5 days past the end to prevent cut-off errors in plots
    tspan = (0.0, total_days + 5.0)

    # 2. Setup Parameters
    if isnothing(fitted_params)
        println("No fitted parameters provided. Using defaults.")
        # Default safe values
        C0, bk, eps = 1.0, 0.005, 800.0
    else
        C0, bk, eps = fitted_params
    end

    params = MosquitoModelDynamics.MosquitoModelParams(
        Nⱼ = Constants.N_TRAPS,
        H₀ = Constants.N_HOUSEHOLDS,
        j  = Constants.ALPHA,
        k  = Constants.K,
        C₀ = C0,
        bₖ = bk,
        ϵ  = eps,
        temp_interp = temp_interp
    )

    # 3. Setup and Solve ODE
    #    Recalculate u0 because it depends on C0
    u0 = MosquitoModelDynamics.default_u0(params)
    prob = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0, tspan, params)

    println("Solving ODE...")
    sol = solve(prob, solver; saveat=saveat_daily, reltol=reltol, abstol=abstol)

    if sol.retcode != :Success
        @warn "Simulation solver failed with status: $(sol.retcode)"
    end

    # 4. Generate Outputs

    # A. Smooth Trajectories (for plotting lines)
    daily_A = [u[MosquitoModelDynamics.IX_A] for u in sol.u]
    daily_M = [u[MosquitoModelDynamics.IX_M] for u in sol.u]

    # Map simulation time t back to real Dates
    sim_dates = [start_date + Day(round(Int, t)) for t in sol.t]

    # B. Discrete Predictions (for scatter plotting against Truth)
    #    We need to calculate model predictions specifically at the rows in trap_df
    data_times = Float64[Dates.value(d - start_date) for d in trap_df.date]

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
        mfai_times = data_times, # The specific times matching trap_df
        mfai_pred = mfai_pred,   # The predictions matching trap_df
        sol_obj = sol,
        params = params
    )
end

end # module