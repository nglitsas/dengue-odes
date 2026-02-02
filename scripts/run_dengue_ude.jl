# scripts/run_dengue_ude_fit.jl

using DengueODES
using DengueODES.DengueUDE
using DataFrames
using CSV
using Plots
using DataInterpolations

println("="^70)
println("DENGUE UDE FITTING - DUAL NEURAL NETWORKS")
println("="^70)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
println("[1] Loading data...")

dengue_df = CSV.read("data/processed/dengue_cases.csv", DataFrame)
temp_df = CSV.read("data/processed/temperature_data.csv", DataFrame)
precip_df = CSV.read("data/processed/precipitation_data.csv", DataFrame)  # NEW

# Create interpolations
temp_interp = LinearInterpolation(temp_df.temperature, temp_df.time)
precip_interp = LinearInterpolation(precip_df.precipitation, precip_df.time)  # NEW

t_start = minimum(dengue_df.time)
observed_times = Float64.(dengue_df.time)
observed_incidence = Float64.(dengue_df.cases)

println("   Observations: $(length(observed_times))")

# ============================================================================
# 2. BUILD DUAL-NN UDE MODEL
# ============================================================================
println("\n[2] Building dual-NN UDE model...")
println("   Physiological NN: Temperature → [βₘ, βₕ, θₘ, b]")
println("   Ecological NN: [Temperature, Precipitation] → C")

container, nn_params_init = build_ude_model(
    temp_interp,
    precip_interp,  # NEW: pass precipitation
    t_start;
    # NN architectures
    physio_hidden_dim = 32,
    physio_hidden_layers = 2,
    eco_hidden_dim = 32,
    eco_hidden_layers = 2,
    seed = 42
)

println("   Total NN parameters: $(length(nn_params_init))")
println("     - Physiological: $(length(nn_params_init.physio))")
println("     - Ecological: $(length(nn_params_init.eco))")

# ============================================================================
# 3. TRAIN
# ============================================================================
u0 = DengueODES.model.default_u0(container.base_params, t_start)
tspan = (observed_times[1], observed_times[end])

println("\n[3] Training...")

nn_params_trained, loss_history = train_ude!(
    container,
    nn_params_init,
    observed_times,
    observed_incidence,
    u0,
    tspan;
    optimizer = Optimisers.ADAM(0.001),
    max_iter = 1000,
    verbose = true,
    save_every = 100
)

# ============================================================================
# 4. ANALYZE RESULTS
# ============================================================================
println("\n[4] Generating predictions and extracting learned functions...")

results = predict_ude(container, nn_params_trained, u0, tspan, observed_times)
learned_rates = extract_learned_rates(container, nn_params_trained, observed_times)

# ============================================================================
# 5. PLOTTING
# ============================================================================
println("\n[5] Creating plots...")

# Plot learned functions split by input type
p1 = plot(title="Physiological Rates (Temperature-driven)")
plot!(p1, observed_times, learned_rates.βₘ, label="βₘ (Mosq→Human)", lw=2)
plot!(p1, observed_times, learned_rates.βₕ, label="βₕ (Human→Mosq)", lw=2)

p2 = plot(title="Physiological Rates (Temperature-driven)")
plot!(p2, observed_times, learned_rates.θₘ, label="θₘ (Extrinsic Incubation)", lw=2, color=:green)
plot!(p2, observed_times, learned_rates.b, label="b (Biting Rate)", lw=2, color=:purple)

p3 = plot(title="Carrying Capacity (Temperature + Precipitation)")
plot!(p3, observed_times, learned_rates.C, label="C (Carrying Capacity)", lw=2, color=:orange)

p4 = plot(observed_times, observed_incidence, label="Observed", marker=:o, alpha=0.6)
plot!(p4, results.times, results.incidence, label="UDE Prediction", lw=2, color=:red)
plot!(p4, title="Model Fit", xlabel="Time", ylabel="Incidence")

p_combined = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 800))
savefig(p_combined, "outputs/disease_ude/dual_nn_results.png")

println("\n" * "="^70)
println("COMPLETE - Dual NN UDE with separate input drivers")
println("="^70)