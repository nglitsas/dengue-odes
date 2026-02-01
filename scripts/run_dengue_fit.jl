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
using DengueODES.Dengue.Data       # get_dengue_data
using DengueODES.Dengue.Fit
using DengueODES.Dengue.Report     # plot_dengue_simulation, summarize_dengue_fit, dengue_incidence_comp_table

mkpath("outputs/dengue_fit")

println("--- Starting Dengue Fit Pipeline ---")

# =======================================================
# 1. LOAD DATA
# =======================================================
println("\n[1] Loading Data...")

# Raw cases path (matches truth script)
raw_cases_path = joinpath(@__DIR__, "..", "data", "raw", "dengue_cases-2010_2022.csv")

# Process cases (uses your real function with mean smoothing)
case_df = Dengue.get_dengue_data(filename=raw_cases_path, mean=true)

fit_start_date = minimum(case_df.dt_sin_pri)
fit_end_date   = maximum(case_df.dt_sin_pri)

weather_df = Temperature.get_weather_data(force_process=true)

# Small buffer for interpolation endpoints
buffer = Day(1)
weather_subset = filter(row -> (fit_start_date - buffer) <= row.date <= (fit_end_date + buffer), weather_df)

temp_interp = Temperature.get_temperature_interpolator(weather_subset)

println("   Case date range: $fit_start_date to $fit_end_date")

# =======================================================
# 2. PREPARE DATA FOR FITTING (GLOBAL TIME AXIS)
# =======================================================
println("\n[2] Preparing vectors for fitting...")

# observed_times in GLOBAL t units (date_to_t)
observed_times = Float64[TimeUtil.date_to_t(d) for d in case_df.dt_sin_pri]

# Incidence column: from your truth script it's :notified
incidence_col = :notified
if hasproperty(case_df, incidence_col)
    observed_incidence = Float64.(case_df[!, incidence_col])
else
    error("case_df missing incidence column: $incidence_col. Check your data.")
end

println("   Fitting window (global t): $(minimum(observed_times)) to $(maximum(observed_times))")

# =======================================================
# 3. RUN PARAMETER FITTING
# =======================================================
println("Julia started with $(Threads.nthreads()) threads available")
println("\n[3] Fitting model parameters...")
training_days = maximum(observed_times) - minimum(observed_times)

fit_result = Fit.fit_dengue_model(
    observed_times,
    observed_incidence,
    temp_interp;                 # GLOBAL-time temp interpolator
    training_days = Float64(training_days),  # use all data for fitting
    config_path = "data/dengue_fit.toml",
    method = :lhs,
    n_lhs_samples = 800,
    solver = Rosenbrock23()
)

fitted_params = fit_result.param
println("   Fitted Params: ", Dict(fit_result.fit_config.names[j] => round(fitted_params[j], sigdigits=6) for j in 1:length(fitted_params)))

if hasproperty(fit_result, :converged)
    # LsqFit result
    println("   Converged:  ", fit_result.converged)
    println("   Residual SSE: ", sum(fit_result.resid .^ 2))
else
    # LHS result
    println("   Residual SSE: ", fit_result.resid)
end

# =======================================================
# 4. RUN SIMULATION & REPORT (WITH FITTED PARAMS)
# =======================================================
println("\n[4] Running simulation & generating report...")

report_results = Report.plot_dengue_simulation(
    case_df,
    temp_interp;
    fitted_params = fitted_params,
    solver = Tsit5()
)

# =======================================================
# 5. SAVE OUTPUTS
# =======================================================
println("\n[5] Saving outputs...")

savefig(report_results.sim_plot, "outputs/dengue_fit/fit_plot.png")
savefig(report_results.resid_plot, "outputs/dengue_fit/fit_resids_plot.png")

stats = Report.summarize_dengue_fit(case_df, report_results.sim)
stats_df = DataFrame(Metric=["SSE","RMSE","MAE"], Value=[stats.SSE, stats.RMSE, stats.MAE])
CSV.write("outputs/dengue_fit/fit_stats.csv", stats_df)

# Save fitted parameters with dynamic names from config
param_names = fit_result.fit_config.names
param_df = DataFrame(Parameter=param_names, Value=fitted_params)
CSV.write("outputs/dengue_fit/fitted_params.csv", param_df)

incidence_df = Report.dengue_incidence_comp_table(case_df, report_results.sim)
CSV.write("outputs/dengue_fit/incidence_observed_vs_pred.csv", incidence_df)

println("   Saved incidence table: outputs/dengue_fit/incidence_observed_vs_pred.csv")
println("   All results saved to outputs/dengue_fit/")
println("\n=== Dengue fit pipeline complete ===")