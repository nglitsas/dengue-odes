# src/models/disease/Fit_Dengue_UDE.jl
module DengueUDE

using DifferentialEquations
using SciMLSensitivity
using Optimization
using OptimizationOptimisers
using OptimizationOptimJL
using ComponentArrays
using Lux
using Random
using Printf
using Statistics
using JLD2
using Dates

using ..model
using ..EpidemEnto

export build_ude_model, train_ude!, predict_ude, save_nn_weights, load_nn_weights

# ============================================================================
# 1. DUAL NEURAL NETWORK ARCHITECTURE
# ============================================================================
"""
Build TWO separate neural networks:
1. Physiological NN: T → [βₘ, βₕ, θₘ, b] (4 outputs)
2. Ecological NN: [T, P] → C (1 output for carrying capacity)
"""
function build_dual_nn(;
                       # Physiological NN (temperature-driven rates)
                       physio_input_dim::Int = 1,    # Just temperature
                       physio_output_dim::Int = 4,   # βₘ, βₕ, θₘ, b
                       physio_hidden_dim::Int = 32,
                       physio_hidden_layers::Int = 2,
                       # Ecological NN (environmental carrying capacity)
                       eco_input_dim::Int = 2,       # Temperature + Precipitation
                       eco_output_dim::Int = 1,      # Carrying capacity multiplier
                       eco_hidden_dim::Int = 32,
                       eco_hidden_layers::Int = 2,
                       activation = tanh)
    
    # Physiological rates network (temperature only)
    physio_layers = []
    push!(physio_layers, Dense(physio_input_dim, physio_hidden_dim, activation))
    for _ in 1:(physio_hidden_layers-1)
        push!(physio_layers, Dense(physio_hidden_dim, physio_hidden_dim, activation))
    end
    push!(physio_layers, Dense(physio_hidden_dim, physio_output_dim))
    physio_nn = Chain(physio_layers...)
    
    # Ecological carrying capacity network (temperature + precipitation)
    eco_layers = []
    push!(eco_layers, Dense(eco_input_dim, eco_hidden_dim, activation))
    for _ in 1:(eco_hidden_layers-1)
        push!(eco_layers, Dense(eco_hidden_dim, eco_hidden_dim, activation))
    end
    push!(eco_layers, Dense(eco_hidden_dim, eco_output_dim))
    eco_nn = Chain(eco_layers...)
    
    return (physio = physio_nn, eco = eco_nn)
end

# ============================================================================
# 2. UDE CONTAINER WITH DUAL NNs
# ============================================================================
struct UDEContainer
    base_params::model.ModelParams
    physio_nn::Any      # NN for physiological rates
    eco_nn::Any         # NN for carrying capacity
    nn_params::ComponentArray  # Combined parameters
    physio_state::Any
    eco_state::Any
    t_start::Float64
    precip_interp::Any  # NEW: precipitation interpolation
    output_scales::NamedTuple
end

# ============================================================================
# 3. UDE DYNAMICS WITH DUAL NNs
# ============================================================================
"""
Modified DengueModel! with dual NNs:
- Physiological NN learns: βₘ, βₕ, θₘ, b from Temperature
- Ecological NN learns: C from Temperature + Precipitation
"""
function UDEModel!(du, u, container::UDEContainer, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u
    p = container.base_params
    
    # Get base entomological rates (development, mortality, oviposition)
    r = EpidemEnto.get_rates(t, p.temp_interp, p)
    
    # Totals
    M = Mₛ + Mₑ + Mᵢ
    H = Hₛ + Hₑ + Hᵢ + Hᵣ
    
    # Get current environmental conditions
    T_current = p.temp_interp(t)
    P_current = container.precip_interp(t)
    
    # --- PHYSIOLOGICAL NN (Temperature → Rates) ---
    T_norm = Float32[(T_current - 25.0) / 10.0]  # Normalize temperature
    
    physio_out, _ = container.physio_nn(T_norm, 
                                        container.nn_params.physio, 
                                        container.physio_state)
    
    # Extract learned physiological rates
    βₘₜ = container.output_scales.βₘ_scale * sigmoid(physio_out[1])
    βₕₜ = container.output_scales.βₕ_scale * sigmoid(physio_out[2])
    θₘₜ = container.output_scales.θₘ_scale * softplus(physio_out[3])
    bₜ = container.output_scales.b_scale * softplus(physio_out[4])
    
    # --- ECOLOGICAL NN (Temperature + Precipitation → Carrying Capacity) ---
    # Normalize inputs
    T_norm_eco = (T_current - 25.0) / 10.0
    P_norm = (P_current - 50.0) / 50.0  # Assuming precip in mm, adjust as needed
    
    eco_input = Float32[T_norm_eco, P_norm]
    
    eco_out, _ = container.eco_nn(eco_input, 
                                   container.nn_params.eco, 
                                   container.eco_state)
    
    # Carrying capacity multiplier (centered around 1.0, can vary ±50%)
    C_multiplier = container.output_scales.C_scale * (1.0 + 0.5 * tanh(eco_out[1]))
    C_val = C_multiplier * EpidemEnto.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)
    
    # --- FORCES OF INFECTION ---
    λₘₕ = bₜ * βₘₜ * (Mᵢ / max(H, 1.0))
    λₕₘ = bₜ * βₕₜ * (Hᵢ / max(H, 1.0))
    
    # --- DYNAMICS ---
    du[model.IX_Mₐ] = p.k * r.δₜ * (1 - Mₐ / C_val) * M - (r.γₘₜ + r.μₐₜ + p.cₐ) * Mₐ
    du[model.IX_Mₛ] = r.γₘₜ * Mₐ - λₕₘ * Mₛ - (r.μₘₜ + p.cₘ) * Mₛ
    du[model.IX_Mₑ] = λₕₘ * Mₛ - (θₘₜ + r.μₘₜ + p.cₘ) * Mₑ
    du[model.IX_Mᵢ] = θₘₜ * Mₑ - (r.μₘₜ + p.cₘ) * Mᵢ
    du[model.IX_Hₛ] = p.μₕ * (H - Hₛ) - λₘₕ * Hₛ
    du[model.IX_Hₑ] = λₘₕ * Hₛ - (p.θₕ + p.μₕ) * Hₑ
    du[model.IX_Hᵢ] = p.θₕ * Hₑ - (p.γₕ + p.μₕ) * Hᵢ
    du[model.IX_Hᵣ] = p.γₕ * Hᵢ - p.μₕ * Hᵣ
    
    return nothing
end

# ============================================================================
# 4. BUILD UDE MODEL WITH DUAL NNs
# ============================================================================
function build_ude_model(temp_interp,
                         precip_interp,  # NEW: precipitation data
                         t_start::Float64;
                         # Base parameters
                         N = 256088.0,
                         μₕ = 3.605e-5,
                         θₕ = 0.125,
                         γₕ = 0.125,
                         k = 0.5,
                         cₐ = 0.0,
                         cₘ = 0.0,
                         C₀ = 0.1061326 * model.N_HOUSEHOLDS,
                         bₖ = 0.453566,
                         ϵ = 4910.463,
                         ϕ = 0.14,
                         # Physiological NN architecture
                         physio_hidden_dim::Int = 32,
                         physio_hidden_layers::Int = 2,
                         # Ecological NN architecture
                         eco_hidden_dim::Int = 32,
                         eco_hidden_layers::Int = 2,
                         # Output scaling
                         βₘ_scale::Float64 = 0.75,
                         βₕ_scale::Float64 = 0.75,
                         θₘ_scale::Float64 = 0.1,
                         b_scale::Float64 = 1.0,
                         C_scale::Float64 = 1.0,
                         seed::Int = 42)
    
    # Create base parameters
    base_params = model.ModelParams(
        N = N, μₕ = μₕ, θₕ = θₕ, γₕ = γₕ,
        k = k, cₐ = cₐ, cₘ = cₘ,
        C₀ = C₀, bₖ = bₖ, ϵ = ϵ, ϕ = ϕ,
        t_start = t_start,
        temp_interp = temp_interp
    )
    
    # Build dual NNs
    rng = Random.default_rng()
    Random.seed!(rng, seed)
    
    nns = build_dual_nn(;
        physio_input_dim = 1,
        physio_output_dim = 4,
        physio_hidden_dim = physio_hidden_dim,
        physio_hidden_layers = physio_hidden_layers,
        eco_input_dim = 2,
        eco_output_dim = 1,
        eco_hidden_dim = eco_hidden_dim,
        eco_hidden_layers = eco_hidden_layers
    )
    
    # Initialize both NNs
    physio_params, physio_state = Lux.setup(rng, nns.physio)
    eco_params, eco_state = Lux.setup(rng, nns.eco)
    
    # Combine parameters into single ComponentArray
    nn_params_combined = ComponentArray(
        physio = ComponentArray(physio_params),
        eco = ComponentArray(eco_params)
    )
    
    # Output scales
    output_scales = (
        βₘ_scale = βₘ_scale,
        βₕ_scale = βₕ_scale,
        θₘ_scale = θₘ_scale,
        b_scale = b_scale,
        C_scale = C_scale
    )
    
    container = UDEContainer(
        base_params,
        nns.physio,
        nns.eco,
        nn_params_combined,
        physio_state,
        eco_state,
        t_start,
        precip_interp,  # Store precipitation interpolation
        output_scales
    )
    
    return container, nn_params_combined
end

# ============================================================================
# 5. LOSS FUNCTION (unchanged logic, uses updated container)
# ============================================================================
function ude_loss(nn_params_new::ComponentArray,
                  container_template::UDEContainer,
                  data_times::AbstractVector,
                  data_incidence::AbstractVector,
                  u0::Vector{Float64},
                  tspan::Tuple{Float64, Float64})
    
    # Update container with new NN params
    container_updated = UDEContainer(
        container_template.base_params,
        container_template.physio_nn,
        container_template.eco_nn,
        nn_params_new,
        container_template.physio_state,
        container_template.eco_state,
        container_template.t_start,
        container_template.precip_interp,
        container_template.output_scales
    )
    
    # Solve ODE
    prob = ODEProblem(UDEModel!, u0, tspan, container_updated)
    
    sol = solve(prob, Tsit5();
                saveat = data_times,
                sensealg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),
                maxiters = 1e6,
                abstol = 1e-6,
                reltol = 1e-6)
    
    if sol.retcode != ReturnCode.Success
        return Inf
    end
    
    # Extract incidence
    pred_incidence = container_updated.base_params.θₕ .* sol[model.IX_Hₑ, :]
    
    # MSE loss with regularization
    loss = mean((pred_incidence .- data_incidence).^2)
    reg = 1e-6 * sum(abs2, nn_params_new)
    
    return loss + reg
end

# ============================================================================
# 6. TRAINING (unchanged)
# ============================================================================
function train_ude!(container::UDEContainer,
                    nn_params_init::ComponentArray,
                    observed_times::Vector{Float64},
                    observed_incidence::Vector{Float64},
                    u0::Vector{Float64},
                    tspan::Tuple{Float64, Float64};
                    optimizer = Optimisers.ADAM(0.001),
                    max_iter::Int = 500,
                    verbose::Bool = true,
                    save_every::Int = 50,
                    output_dir::String = "outputs/disease_ude")
    
    mkpath(output_dir)
    mkpath(joinpath(output_dir, "nn_weights"))
    
    loss_fn = (θ, p) -> ude_loss(θ, container, observed_times, 
                                  observed_incidence, u0, tspan)
    
    iter = 0
    loss_history = Float64[]
    
    callback = function (p, l)
        iter += 1
        push!(loss_history, l)
        
        if verbose && iter % 10 == 0
            @printf("Iter %4d | Loss: %.6e\n", iter, l)
        end
        
        if iter % save_every == 0
            checkpoint_file = joinpath(output_dir, "nn_weights", 
                                      "checkpoint_iter_$(iter).jld2")
            save_nn_weights(p, checkpoint_file)
            if verbose
                println("   Checkpoint saved: $checkpoint_file")
            end
        end
        
        return false
    end
    
    opt_func = OptimizationFunction(loss_fn, Optimization.AutoZygote())
    opt_prob = OptimizationProblem(opt_func, nn_params_init)
    
    println("\n--- Starting UDE Training ---")
    println("   Output directory: $output_dir")
    result = solve(opt_prob, optimizer; callback = callback, maxiters = max_iter)
    
    println("\n--- Training Complete ---")
    @printf("Final Loss: %.6e\n", result.objective)
    
    final_weights_file = joinpath(output_dir, "nn_weights", "final_weights.jld2")
    save_nn_weights(result.u, final_weights_file)
    println("   Final weights saved: $final_weights_file")
    
    save(joinpath(output_dir, "training_history.jld2"),
         "loss_history", loss_history,
         "final_loss", result.objective,
         "iterations", iter)
    
    return result.u, loss_history
end

# ============================================================================
# 7. WEIGHT SAVING/LOADING
# ============================================================================
function save_nn_weights(nn_params::ComponentArray, filepath::String)
    # Save both physiological and ecological NN weights
    save(filepath, 
         "physio_params", Dict(string(k) => v for (k, v) in pairs(nn_params.physio)),
         "eco_params", Dict(string(k) => v for (k, v) in pairs(nn_params.eco)),
         "timestamp", now())
end

function load_nn_weights(filepath::String)
    data = load(filepath)
    
    physio_dict = data["physio_params"]
    eco_dict = data["eco_params"]
    
    return ComponentArray(
        physio = ComponentArray(physio_dict),
        eco = ComponentArray(eco_dict)
    )
end

# ============================================================================
# 8. PREDICTION
# ============================================================================
function predict_ude(container::UDEContainer,
                     nn_params_trained::ComponentArray,
                     u0::Vector{Float64},
                     tspan::Tuple{Float64, Float64},
                     saveat::Vector{Float64})
    
    container_trained = UDEContainer(
        container.base_params,
        container.physio_nn,
        container.eco_nn,
        nn_params_trained,
        container.physio_state,
        container.eco_state,
        container.t_start,
        container.precip_interp,
        container.output_scales
    )
    
    prob = ODEProblem(UDEModel!, u0, tspan, container_trained)
    sol = solve(prob, Tsit5(); saveat = saveat)
    
    pred_incidence = container_trained.base_params.θₕ .* sol[model.IX_Hₑ, :]
    
    return (solution = sol, incidence = pred_incidence, times = saveat)
end

# ============================================================================
# 9. EXTRACT LEARNED FUNCTIONS
# ============================================================================
"""
Extract learned rate functions from trained dual NNs
"""
function extract_learned_rates(container::UDEContainer,
                                nn_params_trained::ComponentArray,
                                times::Vector{Float64})
    
    n = length(times)
    βₘ_vals = zeros(n)
    βₕ_vals = zeros(n)
    θₘ_vals = zeros(n)
    b_vals = zeros(n)
    C_vals = zeros(n)
    
    for (i, t) in enumerate(times)
        T_current = container.base_params.temp_interp(t)
        P_current = container.precip_interp(t)
        
        # Physiological rates (temperature only)
        T_norm = Float32[(T_current - 25.0) / 10.0]
        physio_out, _ = container.physio_nn(T_norm, 
                                            nn_params_trained.physio, 
                                            container.physio_state)
        
        βₘ_vals[i] = container.output_scales.βₘ_scale * sigmoid(physio_out[1])
        βₕ_vals[i] = container.output_scales.βₕ_scale * sigmoid(physio_out[2])
        θₘ_vals[i] = container.output_scales.θₘ_scale * softplus(physio_out[3])
        b_vals[i] = container.output_scales.b_scale * softplus(physio_out[4])
        
        # Carrying capacity (temperature + precipitation)
        T_norm_eco = (T_current - 25.0) / 10.0
        P_norm = (P_current - 50.0) / 50.0
        eco_input = Float32[T_norm_eco, P_norm]
        
        eco_out, _ = container.eco_nn(eco_input, 
                                      nn_params_trained.eco, 
                                      container.eco_state)
        
        C_mult = container.output_scales.C_scale * (1.0 + 0.5 * tanh(eco_out[1]))
        C_vals[i] = C_mult * EpidemEnto.get_carrying_capacity(t, 
                       container.base_params.C₀, container.base_params.bₖ, 
                       container.base_params.ϵ, container.t_start)
    end
    
    return (βₘ = βₘ_vals, βₕ = βₕ_vals, θₘ = θₘ_vals, b = b_vals, C = C_vals)
end

end # module