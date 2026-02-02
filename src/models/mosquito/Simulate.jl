module Simulate

using DataFrames
using OrdinaryDiffEq
using ..Constants
using ..MosquitoModelDynamics
using ...Shared.TimeUtil

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

    println("\n[DEBUG] Inside Simulation Engine")
    if isnothing(fitted_params)
        println("[DEBUG] Status: USING DEFAULTS (fitted_params is nothing)")
    else
        println("[DEBUG] Status: USING FITTED PARAMS: $fitted_params")
    end

    # 1) Params
    if isnothing(fitted_params)
        C0_val  = 2.0
        bk_val  = 0.1
        eps_val = 0.0
    else
        C0_val, bk_val, eps_val = fitted_params
    end

    # 2) Observation times (collection times)
    t_obs = if :t ∈ names(trap_df)
    Float64.(trap_df.t)
    elseif :date ∈ names(trap_df)
        Float64.(TimeUtil.date_to_t.(trap_df.date))
    else
        error("trap_df must have either :t or :date column")
    end

    @assert issorted(t_obs) "observation times must be sorted increasing"



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

    sol = solve(prob, solver; reltol=1e-6, abstol=1e-6, save_everystep=false)

    # 5) Evaluate states at collection times
    # Assumes ordering A, M, T (as in your CaptureModel_Fast)
    A_vals = [sol(t)[1] for t in t_obs]
    M_vals = [sol(t)[2] for t in t_obs]
    T_vals = [sol(t)[3] for t in t_obs]  # cumulative trapped

    # 6) MFAI at collection times: ΔT / N_TRAPS (NO division by dt)
    @assert :n_traps ∈ names(trap_df) "trap_df must include :n_traps after aggregation"
    ntr = Float64.(trap_df.n_traps)

    mfai_pred = zeros(Float64, length(t_obs))
    mfai_pred[1] = 0.0
    for i in 2:length(t_obs)
        ΔT = T_vals[i] - T_vals[i-1]
        mfai_pred[i] = ΔT / ntr[i]   # <-- key change (per-collection)
    end


    return DataFrame(
        t = t_obs,
        mfai_pred = mfai_pred,
        A_population = A_vals,
        M_population = M_vals,
        T_cumulative = T_vals,
    )
end

end # module
