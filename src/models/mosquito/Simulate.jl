module Simulate

using DifferentialEquations
using DataFrames
using DataInterpolations

using ..MosquitoModelDynamics
using ..Constants

export run_mosquito_simulation

function run_mosquito_simulation(trap_df::DataFrame,
                                 temp_interp;
                                 fitted_params::Union{Nothing, Vector{Float64}} = nothing,
                                 initial_adults::Float64 = 0.0,
                                 solver = Tsit5(),
                                 saveat_daily::Float64 = 1.0,
                                 reltol::Float64 = 1e-6,
                                 abstol::Float64 = 1e-6)

    # 1. Parameter Safety Check (THE FIX)
    # If no parameters are provided (fitted_params is nothing), use defaults.
    # Defaults: [C0=5.0, b_cap=0.5, epsilon=909.0 (July 2019)]
    p = isnothing(fitted_params) ? [1.27, 0.3165, 909.0] : fitted_params

    # 2. Time Span
    t_start = minimum(trap_df.t)
    t_end   = maximum(trap_df.t)
    tspan   = (t_start, t_end)

    # 3. Initial Conditions (u0)
    # Heuristic: Estimate hidden states based on observed Adults (A)
    # Ratios: Eggs >> Larvae > Pupae > Adults
    # If initial_adults is 0, we default to a small seed (0.1) to allow growth
    A0 = max(initial_adults, 0.1) 
    
    u0 = [
        A0 * 20.0,  # Eggs (E): Approx 20x adults
        A0 * 10.0,  # Larvae (L): Approx 10x adults
        A0 * 2.0,   # Pupae (P): Approx 2x adults
        A0          # Adults (A): From Data
    ]

    # 4. Define Problem
    # We pass 'p' (which is now guaranteed to be a Vector) to the ODE
    prob = ODEProblem(mosquito_ode!, u0, tspan, (p, temp_interp))

    # 5. Solve
    sol = solve(prob, solver; saveat=saveat_daily, reltol=reltol, abstol=abstol)

    if sol.retcode != :Success && sol.retcode != :Terminated
        println("Warning: Simulation failed with status $(sol.retcode)")
        return nothing
    end

    # 6. Convert to DataFrame
    df_sim = DataFrame(t = sol.t)
    
    # Extract Adult Population (4th state variable)
    df_sim.mfai_pred = [u[4] for u in sol.u] 

    # Add Date column for easier plotting
    # (Optional, but helpful for debugging)
    # df_sim.date = TimeUtil.t_to_date.(df_sim.t) 

    return df_sim
end

end # module