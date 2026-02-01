module Report

using DataFrames
using Plots
using Dates
using DifferentialEquations
using Statistics

# Import siblings
using ..Simulate

export plot_mosquito_simulation, summarize_mosquito_simulation

"""
    plot_mosquito_simulation(trap_df, temp_interp; kwargs...)
"""
function plot_mosquito_simulation(trap_df::DataFrame,
                                  temp_interp; # <--- REMOVED weather_start_date
                                  solver = Rosenbrock23(),
                                  saveat_daily::Float64 = 1.0,
                                  reltol::Float64 = 1e-6,
                                  abstol::Float64 = 1e-6)

    # 1. Run Simulation (No date arg needed anymore)
    sim = Simulate.run_mosquito_simulation(
        trap_df,
        temp_interp;
        solver=solver,
        saveat_daily=saveat_daily,
        reltol=reltol,
        abstol=abstol
    )

    # 2. Extract Data for Plotting
    y_obs = Float64.(trap_df.mfai_obvs)
    y_hat = sim.mfai_pred
    resid = y_obs .- y_hat

    # 3. Create Plots
    # We use the DATE dates for the x-axis so it is readable (e.g., "2017")
    # instead of raw integers like "6000"
    
    # Sim Plot (Scatter data vs Line model)
    sim_plot = scatter(
        trap_df.date, y_obs;
        label="Observed MFAI",
        xlabel="Date",
        ylabel="MFAI",
        title="Mosquito Simulation (MFAI)",
        legend=:topleft,
        alpha=0.6
    )
    # Overlay the smooth simulation line
    # plot!(sim_plot, sim.dates, sim.M; label="Simulated Adult Females (Scaled?)", lw=2)
    
    # Note: If M is population and MFAI is an index, they might have different scales.
    # Usually we compare MFAI_pred vs MFAI_obs.
    # Let's plot the specific predictions at the data points to be precise:
    scatter!(sim_plot, trap_df.date, y_hat; shape=:star, color=:red, label="Predicted MFAI points")

    # Residual Plot
    resid_plot = scatter(
        trap_df.date, resid;
        label="Residual (obs - pred)",
        xlabel="Date",
        ylabel="Residual",
        title="Mosquito Simulation Residuals",
        legend=:topleft
    )
    hline!(resid_plot, [0.0]; label="0", lw=2, color=:black)

    return (sim_plot=sim_plot, resid_plot=resid_plot, sim=sim)
end

"""
    summarize_mosquito_simulation(trap_df, sim)
"""
function summarize_mosquito_simulation(trap_df::DataFrame, sim)
    y_obs = Float64.(trap_df.mfai_obvs)
    y_hat = sim.mfai_pred
    resid = y_obs .- y_hat

    sse  = sum(resid .^ 2)
    rmse = sqrt(mean(resid .^ 2))
    mae  = mean(abs.(resid))

    return (SSE=sse, RMSE=rmse, MAE=mae)
end

end # module



# module Report

# using Plots
# using Statistics
# using DataFrames

# # Import siblings
# using ..Simulate

# export plot_mosquito_fit, summarize_mosquito_fit

# """
#     plot_mosquito_fit(trap_df, temp_interp, fitted_params; kwargs...)

# Generates plots comparing observed MFAI to model-predicted MFAI and residuals.

# Arguments:
# - `trap_df`: DataFrame with `:date` and `:mfai` (or whatever your processed column is)
# - `temp_interp`: temperature interpolation object
# - `fitted_params`: Vector [C₀, bₖ, ϵ]

# Returns:
# - NamedTuple with:
#     - `fit_plot`: plot of observed vs predicted MFAI over time
#     - `resid_plot`: residuals over time
#     - `sim`: simulation results NamedTuple from run_mosquito_simulation
# """
# function plot_mosquito_fit(trap_df::DataFrame,
#                            temp_interp,
#                            fitted_params::Vector{Float64};
#                            solver = Tsit5(),
#                            saveat_daily::Float64 = 1.0,
#                            reltol::Float64 = 1e-5,
#                            abstol::Float64 = 1e-5)

#     # Run forward simulation (single source of truth for predictions)
#     sim = Simulate.run_mosquito_simulation(
#         trap_df,
#         temp_interp,
#         fitted_params;
#         solver=solver,
#         saveat_daily=saveat_daily,
#         reltol=reltol,
#         abstol=abstol
#     )

#     # Observed MFAI values from processed data
#     # CHANGE HERE if your processed column name is different
#     y_obs = Float64.(trap_df.mfai)

#     # Predicted MFAI (aligned with trap_df rows)
#     y_hat = sim.mfai_pred

#     resid = y_obs .- y_hat

#     # --- Plot 1: Observed vs Predicted MFAI ---
#     fit_plot = scatter(
#         sim.mfai_times,
#         y_obs;
#         label="Observed MFAI",
#         xlabel="Days",
#         ylabel="MFAI",
#         title="Mosquito Model Fit (MFAI)",
#         legend=:topleft
#     )

#     plot!(
#         fit_plot,
#         sim.mfai_times,
#         y_hat;
#         label="Predicted MFAI",
#         lw=2
#     )

#     # --- Plot 2: Residuals over Time ---
#     resid_plot = scatter(
#         sim.mfai_times,
#         resid;
#         label="Residual (obs - pred)",
#         xlabel="Days",
#         ylabel="Residual",
#         title="Mosquito Fit Residuals",
#         legend=:topleft
#     )

#     hline!(resid_plot, [0.0]; label="0", lw=2)

#     return (fit_plot=fit_plot, resid_plot=resid_plot, sim=sim)
# end


# """
#     summarize_mosquito_fit(trap_df, sim)

# Computes basic fit statistics using observed MFAI and simulation predictions.
# """
# function summarize_mosquito_fit(trap_df::DataFrame, sim)
#     y_obs = Float64.(trap_df.mfai)
#     y_hat = sim.mfai_pred
#     resid = y_obs .- y_hat

#     sse  = sum(resid .^ 2)
#     rmse = sqrt(mean(resid .^ 2))
#     mae  = mean(abs.(resid))

#     return (
#         SSE = sse,
#         RMSE = rmse,
#         MAE = mae
#     )
# end

# end # module
