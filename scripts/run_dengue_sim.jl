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
using DengueODES.Dengue          # The module defined in DengueTransmission.jl
using DengueODES.Dengue.Report   # The submodule

# Ensure output directory exists
mkpath("outputs/disease_sim")

println("--- Starting Human Dengue Simulation Pipeline ---")

# -------------------------------------------------------
# 1. Load Data (Cases + Weather)
# -------------------------------------------------------

# A. Process Cases (Target Data)
println("Loading dengue case data...")
# Uses the function from data.jl you uploaded
cases_df = DengueODES.Dengue.get_dengue_data() 
# NOTE: Your CSV uses 'dt_sin_pri' for dates, not 'date'
sim_start_date = minimum(cases_df.dt_sin_pri)
sim_end_date   = maximum(cases_df.dt_sin_pri)

println("Simulation Range: $sim_start_date to $sim_end_date")

# B. Process Weather (Driver Data)
println("Loading and interpolating weather...")
weather_df = Temperature.get_weather_data(force_process=false)

# Buffer: Add 7 days before/after to prevent solver edge-cases
buffer = Day(7)
weather_subset = filter(row -> (sim_start_date - buffer) <= row.date <= (sim_end_date + buffer), weather_df)

# Create Interpolator
temp_interp = Temperature.get_temperature_interpolator(weather_subset)

# -------------------------------------------------------
# 2. Run Simulation / Fit
# -------------------------------------------------------
println("Running dengue transmission model...")

# reporting_rate: The fraction of total infections that are clinically reported.
# Start with 0.15 (15%) and adjust if the red line is too high/low.
reporting_rate_guess = 0.15 

results = Report.plot_dengue_simulation(
    cases_df, 
    temp_interp;
    reporting_rate = reporting_rate_guess, 
    solver = Rosenbrock23() # Robust solver for stiff biological models
)

# -------------------------------------------------------
# 3. Save Outputs
# -------------------------------------------------------
println("Saving outputs to outputs/disease_sim/ ...")

# Save the combined plot (Simulation + Residuals)
savefig(results.plot, "outputs/disease_sim/dengue_fit_plot.png")

# Calculate summary stats (RMSE, MAE)
stats = Report.summarize_dengue_fit(cases_df, results.sim, reporting_rate_guess)

# Print stats to console
println("\n--- Fit Statistics ---")
println("RMSE: $(stats.RMSE)")
println("MAE:  $(stats.MAE)")

CSV.write(
    "outputs/disease_sim/sim_stats.csv",
    DataFrame([stats])
)

# -------------------------------------------------------
# 4. Diagnostics: The "Explosion Detector"
# -------------------------------------------------------
println("\n--- Running Diagnostics ---")

# 1. Align Data: Find common dates between Model and Data
common_dates = intersect(cases_df.dt_sin_pri, results.sim.dates)
obs_subset   = filter(row -> row.dt_sin_pri in common_dates, cases_df)

# 2. Extract Predictions for those specific days
sim_indices  = [findfirst(==(d), results.sim.dates) for d in obs_subset.dt_sin_pri]
pred_cases   = results.sim.incidence[sim_indices] .* reporting_rate_guess

# 3. Calculate Residuals
residuals    = abs.(pred_cases .- obs_subset.notified)

# 4. Find the Worst Day
max_val, idx = findmax(residuals)
worst_date   = obs_subset.dt_sin_pri[idx]
worst_obs    = obs_subset.notified[idx]
worst_pred   = pred_cases[idx]

# 5. Check Environment on that day
worst_date_t = DengueODES.Shared.TimeUtil.date_to_t(worst_date)
worst_temp   = temp_interp(worst_date_t)

println("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
println("!!! WORST FIT DETECTED !!!")
println("Date:             $worst_date")
println("Observed Cases:   $worst_obs")
println("Predicted Cases:  $(round(worst_pred, digits=1))")
println("Error (Residual): $(round(max_val, digits=1))")
println("Temp on day:      $(round(worst_temp, digits=1)) °C")
println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

println("Pipeline Complete!")