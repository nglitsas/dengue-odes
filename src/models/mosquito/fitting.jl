module Fitting

using DifferentialEquations
using LsqFit

using ..Constants
using ..MosquitoModelDynamics

export fit_mosquito_model, predict_mfai

# ============================================================================
# 1. HELPER: Update Parameters
# ============================================================================
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
        C₀ = Float64(C0_val),
        bₖ = Float64(bk_val),
        ϵ  = Float64(eps_val),
        temp_interp = base_params.temp_interp
    )
end


# ============================================================================
# 2. HELPER: Solve ODE with New Params
# ============================================================================
"""
    solve_capture(prob_template, p_new; kwargs...)

Remakes the problem with new parameters and a re-calculated initial condition u0.
Returns the solution object or `nothing` if the solver diverged.
"""
function solve_capture(prob_template::ODEProblem,
                       p_new::MosquitoModelDynamics.MosquitoModelParams;
                       solver = Tsit5(),
                       reltol::Float64 = 1e-5, # Tighter tol for fitting
                       abstol::Float64 = 1e-5)
    
    # 1. Recalculate Initial Conditions based on new C0
    #    (If C0 changes, the initial mosquito count A0/M0 should likely scale too)
    u0_new = MosquitoModelDynamics.default_u0(p_new)

    # 2. Remake Problem (Fast, non-allocating)
    prob = remake(prob_template; p=p_new, u0=u0_new)

    # 3. Solve (Verbose=false prevents console spam during optimization)
    sol = solve(prob, solver; reltol=reltol, abstol=abstol, verbose=false)

    # 4. Check for success (The temperature clamp fixes most failures, but just in case)
    return sol.retcode == :Success ? sol : nothing
end

# ============================================================================
# 3. THE MODEL FUNCTION
# ============================================================================
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
                      reltol::Float64 = 1e-5,
                      abstol::Float64 = 1e-5)

    p_new = update_params(base_params, collect(p_fitted))

    sol = solve_capture(prob_template, p_new; solver=solver, reltol=reltol, abstol=abstol)

    # Penalize solver failures so LM rejects these parameters
    if sol === nothing
        return fill(1e12, length(t_steps))
    end

    # IMPORTANT: compute_mfai_theo assumes t_steps are in chronological order
    return MosquitoModelDynamics.compute_mfai_theo(sol, t_steps, trap_idx, N_traps)
end


# ============================================================================
# 4. MAIN FITTING PIPELINE
# ============================================================================
"""
    fit_mosquito_model(x_data, y_data, temp_interp; kwargs...) -> LsqFit.Result

Pipeline entry point. 
Automatically slices data to the first year (training_days) and fits parameters.
Returns the `LsqFit.Result` object containing params and fit statistics.
"""
function fit_mosquito_model(observed_times::Vector{Float64},
                            observed_mfai::Vector{Float64},
                            temp_interp;
                            # --- Fitting Settings ---
                            training_days::Float64 = 365.0,
                            
                            # Initial Guess: [C₀, bₖ, ϵ]
                            # bₖ = 0.0 means "assume constant capacity initially"
                            initial_guess::Vector{Float64} = [1.0, 0.005, 800.0],
                            
                            # --- BOUNDS ---
                            # C₀:   0.1 to 5.0 (Reasonable per-household count)
                            # bₖ:   0.0 to 0.01 (LIMITS GROWTH to ~2.0 per year max)
                            # ϵ:    0 to 1820 (Can start anytime)
                            lower_bounds::Vector{Float64} = [0.1, 0.0, 0.0],
                            upper_bounds::Vector{Float64} = [5.0, 0.01, 2000.0],
                            
                            solver = Tsit5(),
                            reltol::Float64 = 1e-5,
                            abstol::Float64 = 1e-5)

    @assert length(observed_times) == length(observed_mfai) "Time and Data vectors must be same length"
    @assert issorted(observed_times) "Observed times must be sorted"

    # --- 1. SLICE DATA TO TRAINING WINDOW ---
    # We only care about data where time <= training_days (365)
    t0 = observed_times[1]
    train_mask = observed_times .<= (t0 + training_days)
    
    time_train = observed_times[train_mask]
    data_train = observed_mfai[train_mask]
    
    if isempty(time_train)
        error("No data found within the first $training_days days.")
    end

    println("Fitting model on first $(maximum(time_train)) days ($(length(time_train)) points).")

    # --- 2. SETUP PROBLEM TEMPLATE ---
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

    # Only simulate as far as we need to fit (speeds up the process)
    warmup = 30.0
    tspan = (time_train[1] - warmup, maximum(time_train) + 5.0)

    u0_template   = MosquitoModelDynamics.default_u0(base_params)
    prob_template = ODEProblem(MosquitoModelDynamics.CaptureModel!, u0_template, tspan, base_params)

    # --- 3. WRAPPER FOR LSQFIT ---
    # This tells LsqFit: "Given times and params, give me the model prediction"
    model_wrapper(t_steps, p) = predict_mfai(
        t_steps, p, base_params, prob_template;
        trap_idx = MosquitoModelDynamics.IX_T,
        N_traps  = Constants.N_TRAPS,
        solver   = solver,
        reltol   = reltol,
        abstol   = abstol
    )

    # --- 4. RUN OPTIMIZATION ---
    # LsqFit compares 'model_wrapper(time_train)' vs 'data_train'
    fit_result = curve_fit(
        model_wrapper,
        time_train,
        data_train,
        initial_guess;
        lower = lower_bounds,
        upper = upper_bounds
    )

    return fit_result
end

end # module