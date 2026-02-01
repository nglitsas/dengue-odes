# #!/usr/bin/env julia
# using Pkg
# Pkg.activate(@__DIR__ * "/..")

# using DengueODES
# using DengueODES.Models.Mosquito

# using CSV, DataFrames
# using DataInterpolations

# # -------------------------------------------------------
# # Load processed mosquito trap data
# # -------------------------------------------------------

# trap_df = CSV.read(
#     "data/processed/mosq_trapped_model.csv",
#     DataFrame
# )

# # Example expected columns:
# # :day (Float64, days since start)
# # :mfai (Float64, observed MFAI)

# x_data = trap_df.day
# y_data = trap_df.mfai


# # -------------------------------------------------------
# # Build temperature interpolator
# # -------------------------------------------------------

# weather_df = CSV.read(
#     "data/raw/weather-2010_2022.csv",
#     DataFrame
# )

# # Adapt to your actual columns
# temps = weather_df.temperature
# days  = weather_df.day

# temp_interp = LinearInterpolation(temps, days)


# # -------------------------------------------------------
# # Run fit
# # -------------------------------------------------------

# println("Starting mosquito model fit...")

# fit = fit_mosquito_model(
#     x_data,
#     y_data,
#     temp_interp;
#     initial_guess = [30.0, 0.3, 100.0]
# )

# println("Converged: ", fit.converged)
# println("Parameters: ", fit.param)
# println("SSE: ", sum(fit.resid .^ 2))


# # -------------------------------------------------------
# # Save output
# # -------------------------------------------------------

# plots = plot_mosquito_fit(trap_df, temp_interp, fit.param)

# savefig(plots.fit_plot, "outputs/mosquito_fit/fit_plot.png")
# savefig(plots.resid_plot, "outputs/mosquito_fit/residuals_plot.png")

# stats = summarize_mosquito_fit(trap_df, plots.sim)
# CSV.write("outputs/mosquito_fit/fit_summary.csv", DataFrame(stats))
