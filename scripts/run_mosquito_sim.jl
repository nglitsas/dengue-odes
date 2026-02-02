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
mkpath("outputs/mosquito_sim")

println("--- Starting Mosquito Model Simulation Pipeline ---")

# -------------------------------------------------------
# 1. Load Data (Trap + Weather)
# -------------------------------------------------------

# A. Process Weather & Create Interpolator
# This now calls your cleaned Temperature.jl

# B. Process Trap Data
# Assuming Data.get_trap_data() is the correct function name in your current setup
trap_df = DengueODES.MosquitoCapture.Data.get_trap_data()
sim_start_date = minimum(trap_df.date)
sim_end_date   = maximum(trap_df.date)

weather_df = Temperature.get_weather_data(force_process=true)

# C. CUT WEATHER DATA
# We keep 1 day before/after as a safety buffer for the solver
buffer = Day(1)
weather_subset = filter(row -> (sim_start_date - buffer) <= row.date <= (sim_end_date + buffer), weather_df)

# D. Create the Interpolator with the SUBSET
temp_interp = Temperature.get_temperature_interpolator(weather_subset)
println("Interpolator Time Range: ", sim_start_date, " to ", sim_end_date)
println("Value at start: ", temp_interp(DengueODES.Shared.TimeUtil.date_to_t(sim_start_date)))

# using Plots

# # 1. Define the fine-grained time range for plotting
# # We convert the dates to numerical 't' values using your TimeUtil
# t_plot_start = DengueODES.Shared.TimeUtil.date_to_t(sim_start_date)
# t_plot_end   = DengueODES.Shared.TimeUtil.date_to_t(sim_end_date)

# # Create a range of t values (e.g., every 0.1 days for a smooth curve)
# t_range = range(t_plot_start, t_plot_end, length=1000)

# # 2. Evaluate the interpolator at these points
# # This is exactly what the ODE solver "sees" during the simulation
# temp_values = [temp_interp(t) for t in t_range]

# # 3. Convert t_range back to Dates for a readable X-axis
# plot_dates = [DengueODES.Shared.TimeUtil.t_to_date(t) for t in t_range]

# # 4. Create the plot
# p_temp = plot(
#     plot_dates, 
#     temp_values,
#     title  = "Temperature Forcing ($sim_start_date to $sim_end_date)",
#     ylabel = "Temperature (°C)",
#     xlabel = "Date",
#     label  = "Linear Interpolation",
#     lw     = 2,
#     color  = :orange,
#     legend = :outertopright,
#     size   = (900, 400)
# )

# # Optional: Overlay the raw weather station data points to check accuracy
# scatter!(
#     p_temp,
#     weather_subset.date,
#     weather_subset.temp,
#     label = "Weather Station Data",
#     color = :black,
#     markersize = 2,
#     alpha = 0.5
# )

# # Display or save the plot
# display(p_temp)
# savefig("outputs/mosquito_sim/temp_check.png")

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

savefig(results.sim_plot, "outputs/mosquito_sim/sim_plot.png")
savefig(results.resid_plot, "outputs/mosquito_sim/sim_resids_plot.png")
savefig(results.rates_plot, "outputs/mosquito_sim/sim_rates_plot.png")

# Calculate summary stats
stats = Report.summarize_mosquito_fit(trap_df, results.sim)

CSV.write(
    "outputs/mosquito_sim/sim_summary.csv",
    DataFrame([stats])
)

println("Done! Results saved to outputs/mosquito_sim/")


# ... (after savefig lines) ...

# -------------------------------------------------------
# 4. DIAGNOSTIC: Find the Explosion Date
# -------------------------------------------------------
println("\n--- RUNNING DIAGNOSTICS ---")

# 1. Create interpolation from the daily simulation (Length ~2131)
interp_model = LinearInterpolation(results.sim.mfai_pred, results.sim.t)

# 2. Extract model values ONLY at the dates where we have trap data (Length 36)
model_at_traps = [interp_model(t) for t in trap_df.t]

# 3. NOW calculate residuals (Length 36 vs Length 36) -> No DimensionMismatch!
residuals = abs.(model_at_traps .- trap_df.mfai_obvs)

# 4. Find the worst day
max_val, idx = findmax(residuals)
worst_date = trap_df.date[idx]

# 5. Get temperature context
worst_date_t = DengueODES.Shared.TimeUtil.date_to_t(worst_date)
worst_temp = temp_interp(worst_date_t)

println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
println("!!! LARGEST DISCREPANCY DETECTED !!!")
println("Date:             $worst_date")
println("Observed MFAI:    $(round(trap_df.mfai_obvs[idx], digits=4))")
println("Model MFAI:       $(round(model_at_traps[idx], digits=4))")
println("Residual Size:    $(round(max_val, digits=4))")
println("Temp on that day: $(round(worst_temp, digits=2)) °C")
println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")