module Fitting

using DifferentialEquations
using StaticArrays
using LsqFit
using LatinHypercubeSampling
using Base.Threads
using Printf

using ..Constants
using ..MosquitoModelDynamics

export fit_mosquito_model

# ============================================================================
# 1. REFINED PREDICTION
# ============================================================================
function predict_mfai_fast(p_vals, t_steps, t0, temp_interp, prob_template; 
                           solver=Tsit5(), reltol=1e-5, abstol=1e-5)
    
    # 1. Construct Parameters using the Struct
    # We pass C₀ as DENSITY (Mosquitoes per Household).
    # We rely on Entomology.get_carrying_capacity doing the multiplication by N_HOUSEHOLDS.
    p_new = MosquitoModelParams(
        p_vals[1],      # C₀ (Density)
        p_vals[2],      # bₖ
        p_vals[3],      # ϵ
        t0,             # t_start
        temp_interp     # temperature function
    )

    # 2. Re-calculate Initial Conditions (u0)
    # The ODE needs TOTAL mosquitoes, so we must scale Density * N_HOUSEHOLDS here.
    total_capacity = p_vals[1] * Constants.N_HOUSEHOLDS
    
    u0_new = SVector{3}(
        0.85 * total_capacity, # Aquatic: 85% of Capacity
        0.50 * total_capacity, # Adult:   50% of Capacity
        0.0                    # Trapped: Starts at 0
    )

    # 3. Solve
    sol = solve(remake(prob_template; p=p_new, u0=u0_new), 
                solver; 
                saveat=t_steps, 
                reltol=reltol, 
                abstol=abstol, 
                maxiters=1e6)

    if sol.retcode != ReturnCode.Success
        return fill(1e12, length(t_steps))
    end

    # 4. Calculate MFAI using the corrected function from ModelDynamics
    # This handles the (Catch / Traps / dt) calculation correctly.
    return compute_mfai_theo(sol, t_steps)
end

# ============================================================================
# 2. REFINED PARALLEL LHS
# ============================================================================
function fit_lhs(time_train, data_train, t0, temp_interp, prob_template, 
                 lower_bounds, upper_bounds; n_samples=2000)

    n_params = length(lower_bounds)
    plan, _ = LHCoptim(n_samples, n_params, 100) 
    bounds = [(lower_bounds[i], upper_bounds[i]) for i in 1:n_params]
    param_samples = scaleLHC(plan, bounds)
    
    sse_values = fill(1e12, n_samples)

    println("--- Executing LHS on $(Threads.nthreads()) Threads ---")

    @threads for i in 1:n_samples
        p = param_samples[i, :]
        y_hat = predict_mfai_fast(p, time_train, t0, temp_interp, prob_template)
        
        sse = 0.0
        for j in 1:length(data_train)
            diff = data_train[j] - y_hat[j]
            sse += diff * diff
        end
        
        sse_values[i] = isfinite(sse) ? sse : 1e12
    end

    best_idx = argmin(sse_values)
    return (param = param_samples[best_idx, :], resid = sse_values[best_idx])
end

# ============================================================================
# 3. MAIN INTERFACE
# ============================================================================
function fit_mosquito_model(observed_times::Vector{Float64},
                            observed_mfai::Vector{Float64},
                            temp_interp;
                            training_days::Float64 = 365.0,
                            initial_guess = [1.33, 0.31, 900.0],
                            lower_bounds = [0.1, 0.0, 0.0], 
                            upper_bounds = [5.0, 1.5, 1820.0], 
                            method::Symbol = :lhs,
                            n_lhs_samples::Int = 2000)

    # Slice training data
    t0 = observed_times[1]
    train_mask = observed_times .<= (t0 + training_days)
    time_train = observed_times[train_mask]
    data_train = observed_mfai[train_mask]

    # Problem Template
    # We initialize with a safe default. The actual values are overwritten in predict_mfai_fast.
    u0_init = SVector{3}(100.0, 100.0, 0.0) 
    
    tspan = (time_train[1], time_train[end])
    
    # Initialize params with the Struct
    p_init = MosquitoModelParams(initial_guess[1], initial_guess[2], initial_guess[3], t0, temp_interp)

    prob_template = ODEProblem(MosquitoModelDynamics.CaptureModel_Fast, u0_init, tspan, p_init)

    # Define model_wrapper for LsqFit
    model_wrapper(t, p) = predict_mfai_fast(p, t, t0, temp_interp, prob_template)

    if method == :lhs
        best_lhs = fit_lhs(time_train, data_train, t0, temp_interp, prob_template, 
                           lower_bounds, upper_bounds; n_samples=n_lhs_samples)
        
        println("--- Refining LHS result with LsqFit ---")
        result = curve_fit(model_wrapper, time_train, data_train, best_lhs.param; 
                           lower=lower_bounds, upper=upper_bounds)
    else
        result = curve_fit(model_wrapper, time_train, data_train, initial_guess; 
                           lower=lower_bounds, upper=upper_bounds)
    end

    @printf("\nFinal Optimized Parameters:\n")
    @printf("C₀ (Density):     %.3f [Mosquitoes/Household]\n", result.param[1])
    @printf("bₖ (Growth Rate): %.4f\n", result.param[2])
    @printf("ϵ  (Shift Day):   %.1f\n", result.param[3])
    
    return result
end

end # module