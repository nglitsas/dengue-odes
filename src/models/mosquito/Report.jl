module Report

using DataFrames
using Plots
using Dates
using DifferentialEquations
using Statistics
using Printf
using DataInterpolations

using ..Simulate
using ..MosquitoModelDynamics
using ..Constants
# Import TimeUtil so we can convert model 't' back to real dates for the plot
using ...Shared.TimeUtil

export plot_mosquito_simulation, summarize_mosquito_fit, mosquito_mfai_comp_table

function align_simulation_to_observations(sim_df::DataFrame, trap_df::DataFrame)
    # Both now use Global Time (t ~ 6000+)
    interp = LinearInterpolation(sim_df.mfai_pred, sim_df.t)
    obs_t = hasproperty(trap_df, :t) ? trap_df.t : trap_df.time
    model_at_obs = interp.(obs_t)
    return model_at_obs
end

function plot_mosquito_simulation(trap_df::DataFrame,
                                  temp_interp;
                                  fitted_params::Union{Nothing, Vector{Float64}} = nothing,
                                  solver = Tsit5(),
                                  saveat_daily::Float64 = 1.0,
                                  reltol::Float64 = 1e-6,
                                  abstol::Float64 = 1e-6)

    # 1. Run Simulation
    sim = Simulate.run_mosquito_simulation(
        trap_df,
        temp_interp;
        fitted_params = fitted_params,
        solver = solver,
        saveat_daily = saveat_daily,
        reltol = reltol,
        abstol = abstol
    )

    if isnothing(sim)
        return nothing
    end

    # 2. Prepare Data
    y_obs = hasproperty(trap_df, :mfai_obvs) ? Float64.(trap_df.mfai_obvs) : Float64.(trap_df.mfai)
    y_hat_aligned = align_simulation_to_observations(sim, trap_df)

    # Plot 1: Time Series
    sim_plot = scatter(
        trap_df.date, y_obs;
        label="Observed Data",
        xlabel="Date",
        ylabel="MFAI",
        title="Model Fit",
        legend=:topleft,
        color=:blue
    )
    
    # --- PLOTTING FIX ---
    # Since t is ~6000+, we must convert it back to Date objects
    # to overlay correctly with trap_df.date
    sim_dates = TimeUtil.t_to_date.(sim.t)

    plot!(sim_plot, sim_dates, sim.mfai_pred; 
          linewidth=2, color=:red, label="Model (Daily)")

    # Plot 2: Residuals
    residuals_vec = y_obs .- y_hat_aligned
    
    resid_plot = plot(
        trap_df.date, 
        residuals_vec;
        seriestype = :stem,
        label = "Residuals",
        ylabel = "Error",
        color = :purple,
        marker = :circle
    )
    hline!(resid_plot, [0.0]; color=:black, linestyle=:dash, label="")

    return (sim_plot=sim_plot, resid_plot=resid_plot, sim=sim)
end

function summarize_mosquito_fit(trap_df::DataFrame, sim)
    y_obs = hasproperty(trap_df, :mfai_obvs) ? Float64.(trap_df.mfai_obvs) : Float64.(trap_df.mfai)
    y_hat_aligned = align_simulation_to_observations(sim, trap_df)
    
    resid = y_obs .- y_hat_aligned

    sse  = sum(resid .^ 2)
    rmse = sqrt(mean(resid .^ 2))
    mae  = mean(abs.(resid))

    println("\n=== FIT METRICS ===")
    @printf("SSE:  %.4f\n", sse)
    @printf("RMSE: %.4f\n", rmse)
    @printf("MAE:  %.4f\n", mae)

    return (SSE=sse, RMSE=rmse, MAE=mae)
end

function mosquito_mfai_comp_table(trap_df::DataFrame, sim)
    y_obs = hasproperty(trap_df, :mfai_obvs) ? Float64.(trap_df.mfai_obvs) : Float64.(trap_df.mfai)
    y_hat_aligned = align_simulation_to_observations(sim, trap_df)

    df = DataFrame(
        date = trap_df.date,
        mfai_observed = y_obs,
        mfai_predicted = y_hat_aligned,
    )
    df.residual = df.mfai_observed .- df.mfai_predicted
    return df
end

end # module