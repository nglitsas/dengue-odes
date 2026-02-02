module Fit

using DifferentialEquations
using LsqFit
using Printf
using LatinHypercubeSampling
using Base.Threads
using TOML

using ..DengueModel  # Your DengueModel!, ModelParams, default_u0, IX_*
using ..EpidemEnto  # get_rates, get_carrying_capacity

export fit_dengue_model, predict_incidence, compute_sse

# ============================================================================
# 0. CONFIG LOADING (from data/dengue_fit.toml)
# ============================================================================
struct FitConfig
    names::Vector{String}
    initial::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

function load_fit_config(file_path = "data/dengue_fit.toml")
    if !isfile(file_path)
        error("Config file not found: $file_path. Create it with [optimize] section.")
    end
    config = TOML.parsefile(file_path)["optimize"]
    return FitConfig(
        config["names"],
        Float64.(config["initial"]),
        Float64.(config["lower"]),
        Float64.(config["upper"])
    )
end

# ============================================================================
# 1. UPDATE PARAMETERS DYNAMICALLY
# ============================================================================
function update_params(base_params::ModelParams, p_fitted::AbstractVector, fit_names::Vector{String})
    param_dict = Dict{Symbol, Any}(f => getfield(base_params, f) for f in fieldnames(ModelParams))
    
    for (i, name) in enumerate(fit_names)
        sym = Symbol(name)
        if !haskey(param_dict, sym)
            error("Unknown parameter to fit: $name (not in ModelParams)")
        end
        param_dict[sym] = Float64(p_fitted[i])
    end

    return ModelParams(;
        N           = param_dict[:N],
        μₕ          = param_dict[:μₕ],
        θₕ          = param_dict[:θₕ],
        γₕ          = param_dict[:γₕ],
        k           = param_dict[:k],
        cₐ          = param_dict[:cₐ],
        cₘ          = param_dict[:cₘ],
        C₀          = param_dict[:C₀],
        bₖ          = param_dict[:bₖ],
        ϵ           = param_dict[:ϵ],
        ϕ           = param_dict[:ϕ],
        t_start     = base_params.t_start,
        temp_interp = base_params.temp_interp
    )
end

# ============================================================================
# 2. SOLVE ODE HELPER
# ============================================================================
function solve_dengue(prob_template::ODEProblem,
                      p_new::ModelParams;
                      solver = Rosenbrock23(),
                      reltol::Float64 = 1e-7,
                      abstol::Float64 = 1e-7,
                      saveat = nothing)

    u0_new = default_u0(p_new, prob_template.tspan[1])
    prob = remake(prob_template; p = p_new, u0 = u0_new)

    sol = saveat === nothing ?
        solve(prob, solver; reltol=reltol, abstol=abstol, maxiters=10^7, verbose=false) :
        solve(prob, solver; reltol=reltol, abstol=abstol, maxiters=10^7, verbose=false, saveat=saveat)

    if sol.retcode != ReturnCode.Success
        @warn "ODE solve failed" retcode = sol.retcode
        return nothing
    end
    return sol
end

# ============================================================================
# 3. PREDICT INCIDENCE (daily new infections = θₕ * Hₑ)
# ============================================================================
function predict_incidence(t_steps::AbstractVector,
                           p_fitted::AbstractVector,
                           fit_names::Vector{String},
                           base_params::ModelParams,
                           prob_template::ODEProblem;
                           solver = Rosenbrock23(),
                           reltol::Float64 = 1e-7,
                           abstol::Float64 = 1e-7)

    p_new = update_params(base_params, p_fitted, fit_names)
    sol = solve_dengue(prob_template, p_new; solver=solver, reltol=reltol, abstol=abstol, saveat=t_steps)

    if sol === nothing
        return fill(1e12, length(t_steps))
    end

    # New infections per time step = θₕ * Exposed humans
    incidence = p_new.θₕ .* sol[IX_Hₑ, :]  # vectorized

    return incidence
end

# ============================================================================
# 4. LOSS FUNCTION (SSE) — easy to swap later
# ============================================================================
function compute_sse(y_hat::AbstractVector, y_obs::AbstractVector)
    residuals = y_obs .- y_hat
    return sum(residuals .^ 2)
end

function compute_sse(p::AbstractVector,
                     time_train::AbstractVector,
                     data_train::AbstractVector,
                     fit_names::Vector{String},
                     base_params::ModelParams,
                     prob_template::ODEProblem;
                     solver = Rosenbrock23(),
                     reltol::Float64 = 1e-7,
                     abstol::Float64 = 1e-7)

    y_hat = predict_incidence(time_train, p, fit_names, base_params, prob_template;
                              solver=solver, reltol=reltol, abstol=abstol)
    return compute_sse(y_hat, data_train)
end

# ============================================================================
# 5. LHS FITTING (parallel)
# ============================================================================
function fit_lhs(time_train::Vector{Float64},
                 data_train::Vector{Float64},
                 fit_config::FitConfig,
                 base_params::ModelParams,
                 prob_template::ODEProblem;
                 n_samples::Int = 1000,
                 solver = Rosenbrock23(),
                 reltol::Float64 = 1e-7,
                 abstol::Float64 = 1e-7)

    println("\n--- Starting Latin Hypercube Sampling with $(Threads.nthreads()) threads ---")
    println("   Samples: $n_samples")
    println("   Optimizing: ", fit_config.names)

    n_params = length(fit_config.names)
    plan, _ = LHCoptim(n_samples, n_params, 1000)
    bounds = [(fit_config.lower[i], fit_config.upper[i]) for i in 1:n_params]
    param_samples = scaleLHC(plan, bounds)  # n_samples × n_params

    sse_values = Vector{Float64}(undef, n_samples)

    @threads for i in 1:n_samples
        p = param_samples[i, :]
        sse = compute_sse(p, time_train, data_train, fit_config.names, base_params, prob_template;
                          solver=solver, reltol=reltol, abstol=abstol)
        sse_values[i] = isfinite(sse) ? sse : 1e12
    end

    best_idx = argmin(sse_values)
    best_sse = sse_values[best_idx]
    best_params = param_samples[best_idx, :]

    println("\n--- LHS Complete ---")
    @printf("   Best SSE: %.4e\n", best_sse)
    println("   Best params: ", Dict(fit_config.names[j] => round(best_params[j], sigdigits=6) for j in 1:n_params))

    return (param = best_params,
        resid = best_sse,
        all_samples = param_samples,
        all_sse = sse_values,
        fit_config = fit_config)
end

# ============================================================================
# 6. MAIN FITTING FUNCTION
# ============================================================================
function fit_dengue_model(observed_times::Vector{Float64},
                          observed_incidence::Vector{Float64},
                          temp_interp;
                          training_days::Float64 = 365.0,
                          config_path::String = "data/dengue_fit.toml",
                          method::Symbol = :lsqfit,
                          n_lhs_samples::Int = 800,
                          solver = Rosenbrock23(),
                          reltol::Float64 = 1e-7,
                          abstol::Float64 = 1e-7)

    fit_config = load_fit_config(config_path)

    # Slice training data
    t0 = observed_times[1]
    train_mask = observed_times .<= (t0 + training_days)
    time_train = observed_times[train_mask]
    data_train = observed_incidence[train_mask]
    t_start = time_train[1]

    # Base params with temperature driver
    base_params = ModelParams(t_start=t_start, temp_interp = temp_interp)

    # Template ODE problem
    warmup = 0.0
    t_end = maximum(time_train) + 2.0
    tspan = (t_start, t_end + warmup)
    u0_template = default_u0(base_params, t_start)
    prob_template = ODEProblem(DengueModel!, u0_template, tspan, base_params)

    if method == :lhs
        return fit_lhs(time_train, data_train, fit_config, base_params, prob_template;
                       n_samples = n_lhs_samples, solver = solver, reltol = reltol, abstol = abstol)

    elseif method == :lsqfit
        iter_count = 0
        model_wrapper = let
            (t, p) -> begin
                iter_count += 1
                y_hat = predict_incidence(t, p, fit_config.names, base_params, prob_template;
                                          solver = solver, reltol = reltol, abstol = abstol)
                
                if iter_count % 5 == 0
                    @printf("Iter %3d | Params: %s\n", iter_count,
                            Dict(fit_config.names[j] => round(p[j], sigdigits=6) for j in 1:length(p)))
                end
                return y_hat
            end
        end

        println("\n--- Starting LsqFit optimization ---")
        return curve_fit(
            model_wrapper,
            time_train,
            data_train,
            fit_config.initial;
            lower = fit_config.lower,
            upper = fit_config.upper
        )
    else
        error("Unknown method: $method. Use :lhs or :lsqfit")
    end
end

end # module