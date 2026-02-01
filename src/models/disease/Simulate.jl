module Simulate

using DifferentialEquations
using Dates
using DataFrames
using Printf
using Statistics

# Import your Dengue Logic
using ..Dengue # Assuming the parent module is named Dengue
using ...Shared.TimeUtil

export run_dengue_simulation

"""
    run_dengue_simulation(case_df, temp_interp; fitted_params=nothing, kwargs...)

Runs a forward simulation of the human dengue model.
Automatically aligns t=0 with the earliest date in case_df.

Returns a NamedTuple with aligned times, dates, full solution, and predicted incidence.
"""
function run_dengue_simulation(case_df::DataFrame,
                               temp_interp;
                               fitted_params::Union{Nothing, Vector{Float64}} = nothing,
                               solver = Rosenbrock23(), # Stiff solver recommended
                               saveat_daily::Float64 = 1.0,
                               reltol::Float64 = 1e-6,
                               abstol::Float64 = 1e-6)

    println("\n=== Starting Dengue Simulation ===")

    # 1. Setup Time Axis 
    # We use the global TimeUtil to ensure t matches the weather interpolation
    start_date = minimum(case_df.dt_sin_pri)
    end_date   = maximum(case_df.dt_sin_pri)
    
    t_start = TimeUtil.date_to_t(start_date)
    t_end   = TimeUtil.date_to_t(end_date)
    tspan   = (t_start, t_end)

    println("Simulation Range: $start_date to $end_date (t=$t_start to $t_end)")

    # 2. Setup Parameters
    # Initialize default parameters
    p = ModelParams(temp_interp=temp_interp)

    # 2a. Overwrite with fitted params if provided
    # (Example: If you are fitting Reporting Rate or Beta)
    # if !isnothing(fitted_params)
    #     p = ModelParams(..., temp_interp=temp_interp) 
    # end

    # 3. Solve ODE
    u0 = default_u0(p, t_start)
    prob = ODEProblem(DengueModel!, u0, tspan, p)

    println("Solving ODE...")
    sol = solve(prob, solver; saveat=saveat_daily, reltol=reltol, abstol=abstol)

    if sol.retcode != :Success
        @warn "Simulation solver failed with status: $(sol.retcode)"
    end

    # 4. Generate Outputs (Incidence Calculation)
    
    # Map simulation time t back to real Dates
    sim_dates = [TimeUtil.t_to_date(t) for t in sol.t]

    # Calculate "True Incidence" (New Infections per day)
    # In SEIR, New Cases = θₕ * H_e (Flow from Exposed to Infected)
    # H_e is index 6 in your model (M_a, M_s, M_e, M_i, H_s, H_e, H_i, H_r)
    IX_He = 6
    
    # Get H_e for every timestep
    daily_He = [u[IX_He] for u in sol.u]
    
    # Calculate daily new infections: Incidence(t) = θₕ * Hₑ(t)
    # Note: This is "True Infections". Reported cases will be lower (Incidence * ReportingRate)
    daily_incidence = p.θₕ .* daily_He

    return (
        times = sol.t,
        dates = sim_dates,
        incidence = daily_incidence,
        sol = sol,
        p = p # Return params so we know the θ used
    )
end

end # module