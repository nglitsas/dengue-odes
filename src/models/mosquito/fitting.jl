module Fitting

using DifferentialEquations
using LsqFit
using Printf # For pretty printing iterations

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
                       abstol::Float64 = 1e-7)

    u0_new = MosquitoModelDynamics.default_u0(p_new)
    prob = remake(prob_template; p=p_new, u0=u0_new)

    sol = solve(
        prob,
        solver;
        reltol=reltol,
        abstol=abstol,
        maxiters=10^7,
        verbose=false
    )

    if sol.retcode != :Success
        @warn "ODE solve failed" retcode=sol.retcode tspan=prob.tspan u0=u0_new p=p_new
        return nothing
    end
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
    sol = solve_capture(prob_template, p_new; solver=solver, reltol=reltol, abstol=abstol)

    if sol === nothing
        return fill(1e12, length(t_steps))
    end

    return MosquitoModelDynamics.compute_mfai_theo(sol, t_steps, trap_idx, N_traps)
end

# ============================================================================
# 4. MAIN FITTING PIPELINE
# ============================================================================
function fit_mosquito_model(observed_times::Vector{Float64},
                            observed_mfai::Vector{Float64},
                            temp_interp;
                            training_days::Float64 = 365.0,
                            initial_guess::Vector{Float64} = [1.0, 0.3, observed_times[1]+800.0], 
                            lower_bounds::Vector{Float64} = [0.1, 0.0, observed_times[1]],
                            upper_bounds::Vector{Float64} = [5.0, 1.2, observed_times[end] + 2000.0],
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
                      # likely 0.0
    t_end   = maximum(time_train) + 2.0

    tspan = (t_start, t_end + warmup)
    u0_template = MosquitoModelDynamics.default_u0(base_params)
    prob_template = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0_template, tspan, base_params)

    # 3. WRAPPER WITH ITERATION PRINTING
    iter_count = 0
    # The 'let' block keeps 'iter_count' local to this function
    model_wrapper = let
        (t, p) -> begin
            iter_count += 1
            y_hat = predict_mfai(t, p, base_params, prob_template; 
                                 solver=solver, reltol=reltol, abstol=abstol)
            
            # Print every 5th evaluation to avoid screen spam
            if iter_count % 5 == 0
                @printf("Iter %3d | C₀: %.4f | bₖ: %.6f | ϵ: %.2f\n", 
                        iter_count, p[1], p[2], p[3])
            end
            return y_hat
        end
    end

    println("\n--- Starting Optimization ---")
    
  # --- 4. RUN OPTIMIZATION ---
    fit_result = curve_fit(
        model_wrapper,
        time_train,
        data_train,
        initial_guess;
        lower = lower_bounds,
        upper = upper_bounds
        # 'h' is the step size for finite difference Jacobians.
        # Since epsilon (ϵ) is in the hundreds, h=1.0 is appropriate.
        # Note: LsqFit uses the same 'h' for all parameters, 
        # so we'll set it to a small but robust value.
        # h = 1e-4 
    )

    return fit_result
end

end # module