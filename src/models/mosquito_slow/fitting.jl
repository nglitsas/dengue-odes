module Fitting

using DifferentialEquations
using LsqFit
using Printf # For pretty printing iterations 
using LatinHypercubeSampling
using Base.Threads


using ..Constants
using ..MosquitoModelDynamics

export fit_mosquito_model, predict_mfai

# ============================================================================
# 1. HELPER: Update Parameters
# ============================================================================
function update_params(base_params::MosquitoModelDynamics.MosquitoModelParams, p::AbstractVector)
    return MosquitoModelDynamics.MosquitoModelParams(
        Nⱼ = base_params.Nⱼ,
        H₀ = base_params.H₀,
        j  = base_params.j,
        k  = base_params.k,
        C₀ = Float64(p[1]),
        bₖ = Float64(p[2]),
        ϵ  = Float64(p[3]),
        t_start = base_params.t_start,
        temp_interp = base_params.temp_interp
    )
end

# ============================================================================
# 2. HELPER: Solve ODE
# ============================================================================
function solve_capture(prob_template::ODEProblem,
                       p_new::MosquitoModelDynamics.MosquitoModelParams;
                       solver = Rosenbrock23(),
                       reltol::Float64 = 1e-7,
                       abstol::Float64 = 1e-7,
                       saveat = nothing)

    u0_new = MosquitoModelDynamics.default_u0(p_new)
    prob = remake(prob_template; p=p_new, u0=u0_new)

    sol = saveat === nothing ?
        solve(prob, solver; reltol=reltol, abstol=abstol, maxiters=10^7, verbose=false) :
        solve(prob, solver; reltol=reltol, abstol=abstol, maxiters=10^7, verbose=false, saveat=saveat)

    return sol
end

# ============================================================================
# 3. THE MODEL FUNCTION (With Iteration Tracking)
# ============================================================================
function predict_mfai(t_steps::AbstractVector,
                      p_fitted::AbstractVector,
                      base_params::MosquitoModelDynamics.MosquitoModelParams,
                      prob_template::ODEProblem;
                      trap_idx::Int = MosquitoModelDynamics.IX_T,
                      N_traps::Real = Constants.N_TRAPS,
                      solver = Rosenbrock23(),
                      reltol::Float64 = 1e-7,
                      abstol::Float64 = 1e-7)

    p_new = update_params(base_params, p_fitted)
    sol = solve_capture(prob_template, p_new; solver=solver, reltol=reltol, abstol=abstol, saveat=t_steps)
    if sol.retcode != ReturnCode.Success
        return fill(1e12, length(t_steps))
    end

    return MosquitoModelDynamics.compute_mfai_theo(sol, t_steps, trap_idx, N_traps)
end

# ============================================================================
# 4. HELPER: Compute SSE for a parameter set
# ============================================================================
function compute_sse(p::AbstractVector,
                     time_train::AbstractVector,
                     data_train::AbstractVector,
                     base_params::MosquitoModelDynamics.MosquitoModelParams,
                     prob_template::ODEProblem;
                     solver = Rosenbrock23(),
                     reltol::Float64 = 1e-7,
                     abstol::Float64 = 1e-7)
    y_hat = predict_mfai(time_train, p, base_params, prob_template;
                         solver=solver, reltol=reltol, abstol=abstol)
    residuals = data_train .- y_hat
    return sum(residuals .^ 2)
end
# ============================================================================
# 5. LATIN HYPERCUBE SAMPLING METHOD
# ============================================================================
function fit_lhs(time_train::Vector{Float64},
                 data_train::Vector{Float64},
                 base_params::MosquitoModelDynamics.MosquitoModelParams,
                 prob_template::ODEProblem,
                 lower_bounds::Vector{Float64},
                 upper_bounds::Vector{Float64};
                 n_samples::Int = 1000,
                 solver = Rosenbrock23(),
                 reltol::Float64 = 1e-7,
                 abstol::Float64 = 1e-7)

    println("\n--- Starting Latin Hypercube Sampling ---")
    println("   Samples: $n_samples")
    println("   Bounds: C₀ ∈ [$(lower_bounds[1]), $(upper_bounds[1])]")
    println("            bₖ ∈ [$(lower_bounds[2]), $(upper_bounds[2])]")
    println("            ϵ  ∈ [$(lower_bounds[3]), $(upper_bounds[3])]")
    
    # Generate LHS plan
    n_params = length(lower_bounds)
    plan, _ = LHCoptim(n_samples, n_params, 1000)  # 1000 generations for optimization
    
    # Scale samples to parameter bounds
    bounds = [(lower_bounds[i], upper_bounds[i]) for i in 1:n_params]
    param_samples = scaleLHC(plan, bounds)
    
    # ────────────────────────────────────────────────
    # Parallel evaluation
    # ────────────────────────────────────────────────
    sse_values = Vector{Float64}(undef, n_samples)

    @threads for i in 1:n_samples
        p = param_samples[i, :]

        sse = compute_sse(p, time_train, data_train, base_params, prob_template;
                          solver = solver, reltol = reltol, abstol = abstol)

        # Protect against NaN / missing values from failed solves
        sse_values[i] = isfinite(sse) ? sse : 1e12
    end

    # ────────────────────────────────────────────────
    # Find best after all threads are done (serial, very fast)
    # ────────────────────────────────────────────────
    best_idx = argmin(sse_values)
    best_sse = sse_values[best_idx]
    best_params = param_samples[best_idx, :]

    # Optional: collect a few top candidates if you want
    # perm = sortperm(sse_values)
    # top5_idx = perm[1:min(5,length(perm))]
    # ...

    println("\n--- LHS Complete ---")
    @printf("   Best SSE: %.4e\n", best_sse)
    @printf("   Best Params: C₀=%.4f, bₖ=%.6f, ϵ=%.2f\n",
            best_params[1], best_params[2], best_params[3])

    return (param       = best_params,
            resid       = best_sse,
            all_samples = param_samples,
            all_sse     = sse_values)
end
# ============================================================================
# 6. MAIN FITTING PIPELINE (with method selection)  # ← MODIFIED COMMENT
# ============================================================================
function fit_mosquito_model(observed_times::Vector{Float64},
                            observed_mfai::Vector{Float64},
                            temp_interp;
                            training_days::Float64 = 365.0,
                            initial_guess::Vector{Float64} = [1.0, 0.3, observed_times[1]+800.0], 
                            lower_bounds::Vector{Float64} = [0.1, 0.0, observed_times[1]],
                            upper_bounds::Vector{Float64} = [5.0, 1.2, observed_times[end] + 8000.0],
                            method::Symbol = :lsqfit,  # ← ADDED THIS PARAMETER
                            n_lhs_samples::Int = 1200,  # ← ADDED THIS PARAMETER
                            solver = Rosenbrock23(),
                            reltol::Float64 = 1e-7,
                            abstol::Float64 = 1e-7)

    # 1. SLICE DATA
    t0 = observed_times[1]
    train_mask = observed_times .<= (t0 + training_days)
    time_train = observed_times[train_mask]
    data_train = observed_mfai[train_mask]
    t_start = time_train[1]

    # 2. SETUP BASELINE
    base_params = MosquitoModelDynamics.MosquitoModelParams(
        Nⱼ = Constants.N_TRAPS, H₀ = Constants.N_HOUSEHOLDS,
        j  = Constants.ALPHA, k  = Constants.K,
        C₀ = initial_guess[1], bₖ = initial_guess[2], ϵ  = initial_guess[3],
        t_start = t0, 
        temp_interp = temp_interp
    )

    warmup = 0.0
    t_end = maximum(time_train) + 2.0
    tspan = (t_start, t_end + warmup)
    u0_template = MosquitoModelDynamics.default_u0(base_params)
    prob_template = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0_template, tspan, base_params)

    # ← ADDED THIS ENTIRE SECTION (method selection)
    # 3. SELECT METHOD
    if method == :lhs
        # Use Latin Hypercube Sampling
        fit_result = fit_lhs(time_train, data_train, base_params, prob_template,
                            lower_bounds, upper_bounds;
                            n_samples=n_lhs_samples, solver=solver, 
                            reltol=reltol, abstol=abstol)
    
    elseif method == :lsqfit
        # Use LsqFit (original method)
        # ← EVERYTHING BELOW HERE IN THIS BLOCK WAS YOUR ORIGINAL CODE
        iter_count = 0
        model_wrapper = let
            (t, p) -> begin
                iter_count += 1
                y_hat = predict_mfai(t, p, base_params, prob_template; 
                                     solver=solver, reltol=reltol, abstol=abstol)
                
                if iter_count % 5 == 0
                    @printf("Iter %3d | C₀: %.4f | bₖ: %.6f | ϵ: %.2f\n", 
                            iter_count, p[1], p[2], p[3])
                end
                return y_hat
            end
        end

        println("\n--- Starting Optimization (LsqFit) ---")
        
        fit_result = curve_fit(
            model_wrapper,
            time_train,
            data_train,
            initial_guess;
            lower = lower_bounds,
            upper = upper_bounds
        )
    else  # ← ADDED THIS ERROR HANDLING
        error("Unknown fitting method: $method. Use :lsqfit or :lhs")
    end
    # ← END OF ADDED SECTION

    return fit_result
end

end # module