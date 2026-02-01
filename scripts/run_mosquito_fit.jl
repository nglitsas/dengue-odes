#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using CSV
using DataFrames
using DataInterpolations
using Plots
using DifferentialEquations
using Printf

# --- Import Project Modules ---
using DengueODES
using DengueODES.Shared.TimeUtil
using DengueODES.Shared.Temperature  
using DengueODES.MosquitoCapture.Data
using DengueODES.MosquitoCapture.Fitting
using DengueODES.MosquitoCapture.Report 

# Ensure output directory exists
mkpath("outputs/mosquito_fit")

println("--- Starting Mosquito Model Simulation Pipeline ---")

# =======================================================
# 1. LOAD DATA
# =======================================================
println("\n[1] Loading Data...")

# A. Load Weather Data
# We assume this DF has a :date column and a :temp column
weather_df = Temperature.get_weather_data(force_process=false)

# B. Load Trap Data
trap_df = Data.get_trap_data()
start_date = minimum(trap_df.date)
println("   Trap data starts on: $start_date")

# =======================================================
# 2. ALIGN TIME & INTERPOLATORS (CRITICAL STEP)
# =======================================================
println("\n[2] Aligning Time Axes...")

# Base interpolator: global t -> Temp (where t = date_to_t(date))
base_interp = Temperature.get_temperature_interpolator(weather_df)

# Global time for the first trap date, in the SAME units as base_interp
t0 = Float64(date_to_t(start_date))

# Shifted interpolator: local simulation time (days since start_date) -> Temp
shifted_temp_interp(t_sim) = base_interp(t0 + t_sim)

println("   Time alignment complete. t=0.0 corresponds to global t = $t0")

# =======================================================
# 3. PREPARE DATA FOR FITTING
# =======================================================
println("\n[3] Preparing Vectors for Fitting...")

# Convert Trap Dates to Relative Time (Days since start)
observed_times = Float64[Dates.value(d - start_date) for d in trap_df.date]

# Extract counts
if hasproperty(trap_df, :mfai_obvs)
    observed_mfai = Float64.(trap_df.mfai_obvs)
else
    observed_mfai = Float64.(trap_df.mfai)
end

println("   Fitting window: 0.0 to $(maximum(observed_times)) days")

# =======================================================
# 4. RUN PARAMETER FITTING
# =======================================================
println("\n[4] Fitting Model Parameters...")

# We use the SHIFTED interpolator here
fit_result = Fitting.fit_mosquito_model(
    observed_times,
    observed_mfai,
    shifted_temp_interp;  # <--- PASS THE ALIGNED INTERPOLATOR
    training_days = 365.0,
    initial_guess = [1.0, 0.005, 800.0] 
)

fitted_params = fit_result.param
println("   Converged:  ", fit_result.converged)
println("   Fitted Params [C₀, b, ϵ]: ", fitted_params)
println("   Residual SSE: ", sum(fit_result.resid .^ 2))

# =======================================================
# 5. RUN SIMULATION & REPORT
# =======================================================
println("\n[5] Running Simulation & Generating Report...")

# We use the SHIFTED interpolator here too
report_results = Report.plot_mosquito_simulation(
    trap_df, 
    shifted_temp_interp; # <--- PASS THE ALIGNED INTERPOLATOR
    fitted_params = fitted_params, 
    solver = Tsit5()
)

# =======================================================
# 6. SAVE OUTPUTS
# =======================================================
println("\n[6] Saving Outputs...")

# Save Plots
savefig(report_results.sim_plot, "outputs/mosquito_fit/sim_plot.png")
savefig(report_results.resid_plot, "outputs/mosquito_fit/residuals_plot.png")

# Calculate Stats
stats = Report.summarize_mosquito_fit(trap_df, report_results.sim)
stats_df = DataFrame(Metric = ["SSE", "RMSE", "MAE"], Value = [stats.SSE, stats.RMSE, stats.MAE])
CSV.write("outputs/mosquito_fit/sim_stats.csv", stats_df)

# Save Fitted Parameters for record keeping
param_df = DataFrame(Parameter=["C0", "b_cap", "epsilon"], Value=fitted_params)
CSV.write("outputs/mosquito_fit/fitted_params.csv", param_df)

println("   All results saved to outputs/mosquito_fit/")
println("\n=== Pipeline Complete! ===")