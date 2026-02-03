#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using CSV
using DataFrames
using DataInterpolations
using Plots
using DifferentialEquations
using Statistics

# Import your package modules
using DengueODES
using DengueODES.Shared.Temperature 
using DengueODES.Shared.TimeUtil
using DengueODES.Dengue          
using DengueODES.Dengue.Report   
using DengueODES.Dengue.EpidemEnto

# -------------------------------------------------------
# 0. Setup Output Directory
# -------------------------------------------------------
output_dir = joinpath(@__DIR__, "..", "outputs", "dengue_sim")
mkpath(output_dir)
mkpath(joinpath(@__DIR__, "..", "data", "processed"))

println("--- Starting 2019-2020 Dengue Simulation Pipeline ---")
println("Output directory: $output_dir")

# -------------------------------------------------------
# 1. Define Study Period & Load Data
# -------------------------------------------------------
STUDY_START = Date(2019, 6, 1)   # June 1, 2019
STUDY_END = Date(2020, 12, 30)   # December 30, 2020

raw_cases_path = joinpath(@__DIR__, "..", "data", "raw", "dengue_cases-2010_2022.csv")
processed_cases_path = joinpath(@__DIR__, "..", "data", "processed", "dengue_case_data_2019_2020.csv")

println("\n[1] Processing dengue case data...")
println("   Input:  $raw_cases_path")

# Load all data
cases_df_all = Dengue.get_dengue_data(filename=raw_cases_path, mean=true)

# Filter to 2019-2020
cases_df = filter(row -> STUDY_START <= row.dt_sin_pri <= STUDY_END, cases_df_all)

println("   Total cases in full dataset: $(nrow(cases_df_all))")
println("   Cases in 2019-2020 period:   $(nrow(cases_df))")

CSV.write(processed_cases_path, cases_df)

sim_start_date = STUDY_START
sim_end_date = STUDY_END

println("Simulation Range: $sim_start_date to $sim_end_date")

# -------------------------------------------------------
# 2. Load Weather Data
# -------------------------------------------------------
println("\n[2] Loading and interpolating weather...")
weather_df = Temperature.get_weather_data(force_process=false)

# Extended buffer for epsilon parameter (~900 days)
buffer_before = Day(1000)
buffer_after = Day(7)

weather_subset = filter(
    row -> (sim_start_date - buffer_before) <= row.date <= (sim_end_date + buffer_after), 
    weather_df
)

println("   Weather data range: $(minimum(weather_subset.date)) to $(maximum(weather_subset.date))")

temp_interp = Temperature.get_temperature_interpolator(weather_subset)
t0 = TimeUtil.date_to_t(sim_start_date)   # numeric time for June 1, 2019

# Build model params for diagnostics (temp_interp is required field)
p = Dengue.DengueModel.ModelParams(
    temp_interp = temp_interp,
    t_start     = t0        # IMPORTANT: align heaviside / capacity shift to sim window
)

println("Defined t0 = $t0 for start date $sim_start_date")
# Get temperature and rates at start
T_start = temp_interp(t0)
r = EpidemEnto.get_rates(t0, temp_interp, p)
C_val = EpidemEnto.get_carrying_capacity(t0, p.C₀, p.bₖ, p.ϵ, p.t_start)

println("="^60)
println("MOSQUITO VIABILITY TEST")
println("="^60)
println("\nTemperature at t=$t0: $T_start °C")
println("Carrying Capacity: $C_val")

# Test with initial conditions
C_initial = 0.7 * p.N * p.C₀
M_a = 0.35 * C_initial
M_s = 0.50 * C_initial
M_total = M_s  # Assume mostly susceptible

println("\nInitial Mosquitoes:")
println("  Aquatic (Ma): $M_a")
println("  Adults (Ms):  $M_s")
println("  Total:        $(M_a + M_s)")

# Calculate birth vs death
logistic = 1 - (M_a / C_val)
birth_rate = p.k * r.δₜ * logistic * M_total
aquatic_death = (r.γₘₜ + r.μₐₜ + p.cₐ) * M_a

println("\n--- AQUATIC STAGE ---")
println("Oviposition rate (δ):    $(r.δₜ) eggs/female/day")
println("Logistic factor:         $logistic")
println("Total births:            $birth_rate eggs/day")
println("Total deaths+emergence:  $aquatic_death /day")
println("Net change:              $(birth_rate - aquatic_death)")
println("Birth/Death ratio:       $(birth_rate / aquatic_death)")

if birth_rate < aquatic_death
    println("❌ PROBLEM: Births < Deaths! Aquatic stage will crash!")
    println("   Need δ > $(aquatic_death / (p.k * logistic * M_total))")
end

adult_recruitment = r.γₘₜ * M_a
adult_death = (r.μₘₜ + p.cₘ) * M_s

println("\n--- ADULT STAGE ---")
println("Emergence rate (γₘ):     $(r.γₘₜ) /day")
println("Adult mortality (μₘ):    $(r.μₘₜ) /day")
println("Recruitment:             $adult_recruitment adults/day")
println("Deaths:                  $adult_death /day")
println("Net change:              $(adult_recruitment - adult_death)")
println("Recruitment/Death ratio: $(adult_recruitment / adult_death)")

if adult_recruitment < adult_death
    println("❌ PROBLEM: Recruitment < Deaths! Adults will crash!")
end

println("\n--- DIAGNOSIS ---")
if r.δₜ < 1.0
    println("⚠️  Oviposition rate is LOW ($(r.δₜ) eggs/female/day)")
    println("   At $T_start°C, is this expected?")
end

if r.μₘₜ > 0.1
    println("⚠️  Adult mortality is HIGH ($(r.μₘₜ) /day)")
    println("   Lifespan = $(1/r.μₘₜ) days")
end

if r.γₘₜ < 0.05
    println("⚠️  Emergence rate is LOW ($(r.γₘₜ) /day)")
    println("   Development time = $(1/r.γₘₜ) days")
end

println("="^60)

# -------------------------------------------------------
# 3. Run Simulation
# -------------------------------------------------------
println("\n[3] Running dengue transmission model...")

reporting_rate_guess = 0.15 

results = Report.plot_dengue_simulation(
    cases_df, 
    temp_interp;
    reporting_rate = reporting_rate_guess, 
    solver = Rosenbrock23()
)

# -------------------------------------------------------
# 4. Save Main Outputs
# -------------------------------------------------------
println("\n[4] Saving outputs...")

savefig(results.plot, joinpath(output_dir, "dengue_fit_plot.png"))

stats = Report.summarize_dengue_fit(cases_df, results.sim, reporting_rate_guess)

println("\n--- Fit Statistics ---")
println("RMSE: $(stats.RMSE)")
println("MAE:  $(stats.MAE)")

CSV.write(
    joinpath(output_dir, "sim_stats.csv"),
    DataFrame([stats])
)

# -------------------------------------------------------
# 5. Diagnostics
# -------------------------------------------------------
println("\n[5] Running Diagnostics...")

common_dates = intersect(cases_df.dt_sin_pri, results.sim.dates)
obs_subset = filter(row -> row.dt_sin_pri in common_dates, cases_df)

sim_indices = [findfirst(==(d), results.sim.dates) for d in obs_subset.dt_sin_pri]
pred_cases = results.sim.incidence[sim_indices] .* reporting_rate_guess

residuals = abs.(pred_cases .- obs_subset.notified)

max_val, idx = findmax(residuals)
worst_date = obs_subset.dt_sin_pri[idx]
worst_obs = obs_subset.notified[idx]
worst_pred = pred_cases[idx]

worst_date_t = TimeUtil.date_to_t(worst_date)
worst_temp = temp_interp(worst_date_t)

println("\n" * "="^50)
println("WORST FIT DETECTED")
println("="^50)
println("Date:             $worst_date")
println("Observed Cases:   $worst_obs")
println("Predicted Cases:  $(round(worst_pred, digits=1))")
println("Error (Residual): $(round(max_val, digits=1))")
println("Temp on day:      $(round(worst_temp, digits=1))°C")
println("="^50)

# -------------------------------------------------------
# 6. Generate Additional Diagnostic Plots
# -------------------------------------------------------
println("\n[6] Generating diagnostic plots...")

# Temperature overlay
p_temp = Report.plot_temperature_overlay(
    results.sim, cases_df, temp_interp, reporting_rate_guess
)
savefig(p_temp, joinpath(output_dir, "cases_with_temperature.png"))

# Prediction scatter
p_scatter = Report.plot_prediction_scatter(
    results.sim, cases_df, reporting_rate_guess
)
if !isnothing(p_scatter)
    savefig(p_scatter, joinpath(output_dir, "predicted_vs_observed_scatter.png"))
end

# Monthly comparison
p_monthly = Report.plot_monthly_comparison(
    results.sim, cases_df, reporting_rate_guess
)
savefig(p_monthly, joinpath(output_dir, "monthly_comparison.png"))

# Compartment dynamics
p_human = Report.plot_human_compartments(results.sim)
if !isnothing(p_human)
    savefig(p_human, joinpath(output_dir, "human_compartments.png"))
end

p_mosq = Report.plot_mosquito_compartments(results.sim)
if !isnothing(p_mosq)
    savefig(p_mosq, joinpath(output_dir, "mosquito_compartments.png"))
end

println("\n" * "="^50)
println("PIPELINE COMPLETE")
println("="^50)
println("All outputs saved to: $output_dir")
println("="^50)