module Fitting

using DifferentialEquations
using LsqFit

using ..Constants
using ..MosquitoModelDynamics

export fit_mosquito_model, predict_mfai

"""
    update_params(base_params, p)

Returns a new MosquitoModelParams with fixed fields copied from base_params
and fitted fields replaced with p = [C₀, bₖ, ϵ].
"""
function update_params(base_params::MosquitoModelDynamics.MosquitoModelParams, p::AbstractVector)
    @assert length(p) == 3 "Expected p = [C₀, bₖ, ϵ]"
    C0_val, bk_val, eps_val = p

    return MosquitoModelDynamics.MosquitoModelParams(
        Nⱼ = base_params.Nⱼ,
        H₀ = base_params.H₀,
        j  = base_params.j,
        k  = base_params.k,
        C₀ = C0_val,
        bₖ = bk_val,
        ϵ  = eps_val,
        temp_interp = base_params.temp_interp
    )
end


"""
    solve_capture(prob_template, p_new; kwargs...) -> sol or nothing

Remakes prob_template with new parameters AND new u0 (since u0 depends on C₀),
solves, and returns sol if successful else nothing.
"""
function solve_capture(prob_template::ODEProblem,
                       p_new::MosquitoModelDynamics.MosquitoModelParams;
                       solver = Tsit5(),
                       reltol::Float64 = 1e-4,
                       abstol::Float64 = 1e-4)

    u0_new = MosquitoModelDynamics.default_u0(p_new)
    prob   = remake(prob_template; p=p_new, u0=u0_new)

    sol = solve(prob, solver; reltol=reltol, abstol=abstol, verbose=false)
    return sol.retcode == :Success ? sol : nothing
end


"""
    predict_mfai(t_steps, p_fitted, base_params, prob_template; kwargs...) -> yhat

Returns predicted MFAI values at t_steps for parameter vector p_fitted.

This is the function you pass to `LsqFit.curve_fit` as model(x, p).
Uses MosquitoModelDynamics.compute_mfai_theo.
"""
function predict_mfai(t_steps::AbstractVector,
                      p_fitted::AbstractVector,
                      base_params::MosquitoModelDynamics.MosquitoModelParams,
                      prob_template::ODEProblem;
                      trap_idx::Int = MosquitoModelDynamics.IX_T,
                      N_traps::Real = Constants.N_TRAPS,
                      solver = Tsit5(),
                      reltol::Float64 = 1e-4,
                      abstol::Float64 = 1e-4)

    p_new = update_params(base_params, collect(p_fitted))

    sol = solve_capture(prob_template, p_new; solver=solver, reltol=reltol, abstol=abstol)

    # Penalize solver failures so LM rejects these parameters
    if sol === nothing
        return fill(Inf, length(t_steps))
    end

    # IMPORTANT: compute_mfai_theo assumes t_steps are in chronological order
    return MosquitoModelDynamics.compute_mfai_theo(sol, t_steps, trap_idx, N_traps)
end


"""
    fit_mosquito_model(x_data, y_data, temp_interp; kwargs...) -> LsqFit.Result

Pipeline entry point. Assumes x_data/y_data were prepared elsewhere.

Calling this function WILL run the Levenberg–Marquardt fit.
Just importing/including this module will NOT run it.
"""
function fit_mosquito_model(x_data::Vector{Float64},
                            y_data::Vector{Float64},
                            temp_interp;
                            initial_guess::Vector{Float64} = [30.0, 0.3, 100.0],
                            lower_bounds::Vector{Float64} = [1.0, 0.001, 0.0],
                            upper_bounds::Vector{Float64} = [1000.0, 10.0, 2000.0],
                            solver = Tsit5(),
                            reltol::Float64 = 1e-4,
                            abstol::Float64 = 1e-4)

    @assert length(x_data) == length(y_data) "x_data and y_data must have same length"
    @assert issorted(x_data) "x_data must be sorted increasing (MFAI uses deltas between consecutive times)"

    # ---------------------------
    # DATA PIPELINE STUB (external)
    # ---------------------------
    # x_data, y_data = load_processed_trap_series(...)
    # temp_interp    = build_temperature_interpolator(...)
    # ---------------------------

    base_params = MosquitoModelDynamics.MosquitoModelParams(
        Nⱼ = Constants.N_TRAPS,
        H₀ = Constants.N_HOUSEHOLDS,
        j  = Constants.ALPHA,
        k  = Constants.K,
        C₀ = initial_guess[1],
        bₖ = initial_guess[2],
        ϵ  = initial_guess[3],
        temp_interp = temp_interp
    )

    tspan = (0.0, maximum(x_data) + 5.0)

    u0_template   = MosquitoModelDynamics.default_u0(base_params)
    prob_template = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0_template, tspan, base_params)

    model_wrapper(t_steps, p) = predict_mfai(
        t_steps, p, base_params, prob_template;
        trap_idx = MosquitoModelDynamics.IX_T,
        N_traps  = Constants.N_TRAPS,
        solver   = solver,
        reltol   = reltol,
        abstol   = abstol
    )

    # This call runs LM fitting
    fit_result = curve_fit(
        model_wrapper,
        x_data,
        y_data,
        initial_guess;
        lower = lower_bounds,
        upper = upper_bounds
    )

    return fit_result
end

end # module
