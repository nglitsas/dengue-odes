module Simulate

using DifferentialEquations
using DifferentialEquations: ReturnCode  # Add this import
using Dates
using DataFrames
using Printf
using Statistics

# Import your Dengue Logic
using ..Dengue
using ...Shared.TimeUtil
using ..DengueModel
using ..DengueModel: ModelParams, default_u0, DengueModel!
using ..DengueModel: IX_Mₐ, IX_Mₛ, IX_Mₑ, IX_Mᵢ, IX_Hₛ, IX_Hₑ, IX_Hᵢ, IX_Hᵣ

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
                               solver = Rosenbrock23(),
                               saveat_daily::Float64 = 1.0,
                               reltol::Float64 = 1e-8,  # Tighter tolerance
                               abstol::Float64 = 1e-8)  # Tighter tolerance

    println("\n=== Starting Dengue Simulation ===")

    # 1. Setup Time Axis 
    start_date = minimum(case_df.dt_sin_pri)
    end_date   = maximum(case_df.dt_sin_pri)
    
    t_start = TimeUtil.date_to_t(start_date)
    t_end   = TimeUtil.date_to_t(end_date)
    tspan   = (t_start, t_end)

    println("Simulation Range: $start_date to $end_date (t=$t_start to $t_end)")

    # 2. Setup Parameters
    p = ModelParams(temp_interp=temp_interp)

    # 2a. Overwrite with fitted params if provided
    # if !isnothing(fitted_params)
    #     p = ModelParams(..., temp_interp=temp_interp) 
    # end

    # 3. Initial Conditions
    u0 = default_u0(p, t_start)
    
    # 4. Create and Solve ODE Problem
    prob = ODEProblem(DengueModel!, u0, tspan, p)

    println("Solving ODE...")
    sol = solve(
        prob, 
        solver; 
        saveat = saveat_daily, 
        reltol = reltol, 
        abstol = abstol,
        maxiters = 1e7,
        isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u)  # Stop if negative
    )

    # 5. Check Solution Status
    if sol.retcode != ReturnCode.Success
        @warn "Simulation solver failed with status: $(sol.retcode)"
        # Return empty/safe values
        return (
            dates = Date[],
            times = Float64[],
            incidence = Float64[],
            Hs = Float64[], 
            He = Float64[], 
            Hi = Float64[], 
            Hr = Float64[],
            Ma = Float64[], 
            Ms = Float64[], 
            Me = Float64[], 
            Mi = Float64[],
            sol = nothing,
            p = p
        )
    end

    # 6. Extract All Compartments
    Ma_vals = [u[IX_Mₐ] for u in sol.u]
    Ms_vals = [u[IX_Mₛ] for u in sol.u]
    Me_vals = [u[IX_Mₑ] for u in sol.u]
    Mi_vals = [u[IX_Mᵢ] for u in sol.u]

    Hs_vals = [u[IX_Hₛ] for u in sol.u]
    He_vals = [u[IX_Hₑ] for u in sol.u]
    Hi_vals = [u[IX_Hᵢ] for u in sol.u]
    Hr_vals = [u[IX_Hᵣ] for u in sol.u]

    # 7. Calculate Incidence (new infections per day)
    # Incidence = θₕ * Hₑ (flow from Exposed to Infected)
    incidence = p.θₕ .* He_vals

    # 8. Map simulation times back to dates
    sim_dates = [TimeUtil.t_to_date(t) for t in sol.t]

    # 9. Diagnostic: Check for negative values
    if any(x -> x < 0, Ma_vals) || any(x -> x < 0, Hs_vals)
        @warn "Negative populations detected in solution"
        neg_indices = findall(x -> any(val -> val < 0, x), sol.u)
        if !isempty(neg_indices)
            first_neg = neg_indices[1]
            @warn "First negative at t=$(sol.t[first_neg]), date=$(sim_dates[first_neg])"
        end
    end

    println("Simulation complete. Generated $(length(sim_dates)) daily values.")

    return (
        times = sol.t,
        dates = sim_dates,
        incidence = incidence,
        Hs = Hs_vals,
        He = He_vals,
        Hi = Hi_vals,
        Hr = Hr_vals,
        Ma = Ma_vals,
        Ms = Ms_vals,
        Me = Me_vals,
        Mi = Mi_vals,
        sol = sol,
        p = p
    )
end

end # module