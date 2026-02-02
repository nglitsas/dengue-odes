#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using CSV
using DataFrames
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

mkpath("outputs/mosquito_fit")

println("--- Starting Mosquito Fit Pipeline ---")

# =======================================================
# 1. LOAD DATA
# =======================================================
println("\n[1] Loading Data...")

trap_df = Data.get_trap_data()
fit_start_date = minimum(trap_df.date)
fit_end_date   = maximum(trap_df.date)

weather_df = Temperature.get_weather_data(force_process=true)

# --- FIX: INCREASE BUFFER ---
# We need a larger buffer (e.g., 1 year) so the model doesn't crash 
# if the optimizer tests a lag (epsilon) or if we add a warmup period later.
buffer = Year(2)  # Increase from 1 year to 2 years 
weather_subset = filter(row -> (fit_start_date - buffer) <= row.date <= (fit_end_date + buffer), weather_df)

temp_interp = Temperature.get_temperature_interpolator(weather_subset)

println("   Trap date range: $fit_start_date to $fit_end_date")

# =======================================================
# 2. PREPARE DATA
# =======================================================
println("\n[2] Preparing vectors for fitting...")

observed_times = Float64[TimeUtil.date_to_t(d) for d in trap_df.date]

if hasproperty(trap_df, :mfai_obvs)
    observed_mfai = Float64.(trap_df.mfai_obvs)
elseif hasproperty(trap_df, :mfai)
    observed_mfai = Float64.(trap_df.mfai)
else
    error("trap_df missing observed MFAI column (:mfai_obvs or :mfai)")
end

println("   Fitting window (global t): $(minimum(observed_times)) to $(maximum(observed_times))")

# =======================================================
# 3. RUN PARAMETER FITTING
# =======================================================
println("Julia started with $(Threads.nthreads()) threads available")
println("\n[3] Fitting model parameters...")
training_days = maximum(observed_times) - minimum(observed_times)

# --- FIX: LOGICAL BOUNDS FOR EPSILON ---
# Epsilon is a lag/shift (days). It should be small (0-100), not 7000+.
fit_result = Fitting.fit_mosquito_model(
    observed_times,
    observed_mfai,
    temp_interp;                 
    training_days = Float64(training_days),
    # Now C₀ is in units of TOTAL mosquitoes (10^5)
    initial_guess = [1.5e5, 0.3, 900.0],
    lower_bounds  = [1e5, 0.1, 700.0],
    upper_bounds  = [4e5, 0.5, 1100.0],  # ← Increase C₀ max from 2.5e5 to 4e5
    method = :lhs, 
)

fitted_params = fit_result.param
println("   Fitted Params [C₀, b_cap, ϵ]: ", fitted_params)

if hasproperty(fit_result, :converged)
    println("   Converged:   ", fit_result.converged)
    println("   Residual SSE: ", sum(fit_result.resid .^ 2))
else
    println("   Residual SSE: ", fit_result.resid)
end

# =======================================================
# 4. RUN SIMULATION & REPORT
# =======================================================
println("\n[4] Running simulation & generating report...")

report_results = Report.plot_mosquito_simulation(
    trap_df, 
    temp_interp,
    # --- FIX: CORRECT VARIABLE NAME ---
    fitted_params = fitted_params, 
    solver = Tsit5()
)

# =======================================================
# 5. SAVE OUTPUTS
# =======================================================
println("\n[5] Saving outputs...")

savefig(report_results.sim_plot, "outputs/mosquito_fit/fit_plot.png")
savefig(report_results.resid_plot, "outputs/mosquito_fit/fit_resids_plot.png")

stats = Report.summarize_mosquito_fit(trap_df, report_results.sim)
# Ensure SSE is accessed correctly (case sensitivity check)
stats_df = DataFrame(Metric=["SSE","RMSE","MAE"], Value=[stats.SSE, stats.RMSE, stats.MAE])
CSV.write("outputs/mosquito_fit/fit_stats.csv", stats_df)

param_df = DataFrame(Parameter=["C0","b_cap","epsilon"], Value=fitted_params)
CSV.write("outputs/mosquito_fit/fitted_params.csv", param_df)

mfai_df = Report.mosquito_mfai_comp_table(trap_df, report_results.sim)
CSV.write("outputs/mosquito_fit/mfai_observed_vs_pred.csv", mfai_df)

println("   Saved MFIA table: outputs/mosquito_fit/mfai_observed_vs_pred.csv")
println("   All results saved to outputs/mosquito_fit/")
println("\n=== Fit pipeline complete ===")
