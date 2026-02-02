module Simulate

using DataFrames
using OrdinaryDiffEq
using ..Constants
using ..MosquitoModelDynamics
using ...Shared.TimeUtil
using ..Entomology

export run_mosquito_simulation

"""
Runs the mosquito capture model and returns predictions at *collection times*.

Key idea:
- Solve once over [t_start, t_end]
- Evaluate cumulative trapped T(t) at the observation times trap_df.t
- Compute MFAI via: MFAI_i = (T(t_i) - T(t_{i-1})) / N_TRAPS
"""
function run_mosquito_simulation(trap_df, temp_interp;
                                 fitted_params=nothing,
                                 solver=Tsit5())

    # println("\n[DEBUG] Inside Simulation Engine")

    # 2) Observation times (collection times)
    # t_obs = trap_df.date
    t_obs = if "t" ∈ names(trap_df)
    Float64.(trap_df.t)
    elseif "date" ∈ names(trap_df)
        Float64.(TimeUtil.date_to_t.(trap_df.date))
    else
        error("trap_df must have either :t or :date column")
    end

    @assert issorted(t_obs) "observation times must be sorted increasing"
    t_start = minimum(t_obs)
    t_span = (t_start, maximum(t_obs))

        # 1) Params
    if isnothing(fitted_params)
        # SHOULD ALL BE RELATIVE PER HOUSEHOLD VALUES
        C0_val  = 1.33
        bk_val  = 0.03165
        eps_val = 909.0
        println("[DEBUG] Status: Fitted params empty. Using DEFAULTS:")
        println("        C0:  $C0_val")
        println("        bk:  $bk_val")
        println("        eps: $eps_val")
    else
        C0_val, bk_val, eps_val = fitted_params
        println("[DEBUG] Status: USING FITTED PARAMS: $fitted_params")
    end

    # 3) Build parameter struct
    p = MosquitoModelParams(
        C0_val,
        bk_val,
        eps_val,
        t_start,
        temp_interp
    )

    # 4) Solve ODE (your model is out-of-place: returns SVector)
    u0 = MosquitoModelDynamics.default_u0(p)
    prob = ODEProblem(MosquitoModelDynamics.CaptureModel, u0, t_span, p)

    sol = solve(prob, solver; saveat=1.0, reltol=1e-6, abstol=1e-6, save_everystep=false)
    # After creating sol, before making df_rates:
    println("\n" * "="^80)
    println("SOLVER TIME POINTS DIAGNOSTIC")
    println("="^80)
    println("Number of time points: $(length(sol.t))")
    println("Time span: $(minimum(sol.t)) to $(maximum(sol.t))")
    println("First 10 time points: $(sol.t[1:min(10, length(sol.t))])")
    println("Time steps (Δt): $(diff(sol.t[1:min(11, length(sol.t))]))")
    println("\nPoints per year: $(length(sol.t) / ((maximum(sol.t) - minimum(sol.t))/365))")
    println("="^80)

    # 5) Generate Detailed Rates DataFrame for all time steps
    # We loop through sol.t (every step the solver took or saved)
    rates_data = map(sol.t) do t
        δ, γ, μ_a, μ_m = Entomology.get_rates(t, p.temp_interp)
        C = Entomology.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)
        return (t=t, oviposition_rate=δ, emergence_rate=γ, 
                mu_a=μ_a, mu_m=μ_m, capacity=C)
    end
    df_rates = DataFrame(rates_data)
    obs_times = [6210.0, 6269.0, 6330.0, 7000.0, 7114.0, 7200.0, 8330.0]
    for t in obs_times
        C = get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)
        println("t=$t → C=$C")
    end

    # 5) Evaluate states at collection times
    # Assumes ordering A, M, T (as in your CaptureModel_Fast)
    A_vals = [sol(t)[1] for t in t_obs]
    M_vals = [sol(t)[2] for t in t_obs]
    T_vals = [sol(t)[3] for t in t_obs]  # cumulative trapped

    # 6) MFAI at collection times: ΔT / N_TRAPS (NO division by dt)
    mfai_pred = MosquitoModelDynamics.compute_mfai_theo(sol, t_obs)


    df_pred = DataFrame(
        t = t_obs,
        df_rates = df_rates,
        mfai_pred = mfai_pred,
        A_population = A_vals,
        M_population = M_vals,
        T_cumulative = T_vals,
    )
    return (
        sim = df_pred,
        rates = df_rates,
    )
end

end # module
