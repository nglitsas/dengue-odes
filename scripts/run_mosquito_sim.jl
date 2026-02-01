#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using CSV
using DataFrames
using DataInterpolations
using Plots
using DifferentialEquations

# Import your package modules
using DengueODES
using DengueODES.Shared.Temperature  # This works now!
using DengueODES.Shared.TimeUtil
using DengueODES.MosquitoCapture
using DengueODES.MosquitoCapture.Report 

# Ensure output directory exists
mkpath("outputs/mosquito_fit")

println("--- Starting Mosquito Model Simulation Pipeline ---")

# -------------------------------------------------------
# 1. Load Data (Trap + Weather)
# -------------------------------------------------------

# A. Process Weather & Create Interpolator
# This now calls your cleaned Temperature.jl
weather_df = Temperature.get_weather_data(force_process=true)
# NOTE: We don't need the second return value (start_date) anymore!
temp_interp = Temperature.get_temperature_interpolator(weather_df)

# B. Process Trap Data
# Assuming Data.get_trap_data() is the correct function name in your current setup
trap_df = DengueODES.MosquitoCapture.Data.get_trap_data()

# -------------------------------------------------------
# 2. Run Simulation / Fit
# -------------------------------------------------------
println("Running mosquito simulation...")

# Note: Ensure your plot_mosquito_simulation function accepts the interpolator
results = Report.plot_mosquito_simulation(
    trap_df, 
    temp_interp, 
    solver = Rosenbrock23() # <--- CHANGED from default Tsit5()
)

# -------------------------------------------------------
# 3. Save Outputs
# -------------------------------------------------------
println("Saving outputs...")

savefig(results.sim_plot, "outputs/mosquito_fit/sim_plot.png")
savefig(results.resid_plot, "outputs/mosquito_fit/residuals_plot.png")

# Calculate summary stats
stats = summarize_mosquito_simulation(trap_df, results.sim)

CSV.write(
    "outputs/mosquito_fit/sim_summary.csv",
    DataFrame([stats])
)

println("Done! Results saved to outputs/mosquito_fit/")


# ... (after savefig lines) ...

# === DIAGNOSTIC: Find the Explosion Date ===
# FIX: Access .sim.mfai_pred (nested inside results)
residuals = abs.(results.sim.mfai_pred .- trap_df.mfai_obvs)

max_val, idx = findmax(residuals)
worst_date = trap_df.date[idx]

# Get the temperature for that specific day
worst_date_t = DengueODES.Shared.TimeUtil.date_to_t(worst_date)
worst_temp = temp_interp(worst_date_t)

println("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
println("!!! EXPLOSION DETECTED !!!")
println("Worst Date:       $worst_date")
println("Residual Size:    $max_val")
println("Temp on that day: $worst_temp °C")
println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")