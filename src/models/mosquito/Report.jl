module Report

using DataFrames
using Plots
using Dates
using DifferentialEquations
using Statistics
using Printf

# Import siblings
using ..Simulate
using ..MosquitoModelDynamics
using ..Constants

export plot_mosquito_simulation, summarize_mosquito_fit, mosquito_mfai_comp_table


"""
    plot_mosquito_simulation(trap_df, temp_interp; fitted_params=nothing, kwargs...)

Runs the mosquito simulation using the provided parameters and generates diagnostic plots.
"""
function plot_mosquito_simulation(trap_df::DataFrame,
                                  temp_interp;
                                  fitted_params::Union{Nothing, Vector{Float64}} = nothing, # <--- ADDED THIS
                                  solver = Tsit5(),
                                  saveat_daily::Float64 = 1.0,
                                  reltol::Float64 = 1e-6,
                                  abstol::Float64 = 1e-6)

    # 1. Run Simulation
    # We must pass 'fitted_params' down to the simulation!
    sim = Simulate.run_mosquito_simulation(
        trap_df,
        temp_interp;
        fitted_params = fitted_params, # <--- IMPORTANT: Pass it here
        solver = solver,
        saveat_daily = saveat_daily,
        reltol = reltol,
        abstol = abstol
    )

    # 2. Extract Data for Plotting
    # Handle column name variations (:mfai vs :mfai_obvs)
    if hasproperty(trap_df, :mfai_obvs)
        y_obs = Float64.(trap_df.mfai_obvs)
    else
        y_obs = Float64.(trap_df.mfai)
    end
    
    y_hat = sim.mfai_pred
    resid = y_obs .- y_hat

    # =========================================================
    # PLOT 1: OBSERVED vs PREDICTED
    # =========================================================
    
    sim_plot = scatter(
        trap_df.date, y_obs;
        label="Observed Data",
        xlabel="Date",
        ylabel="MFAI (Mosquitoes/Trap)",
        title="Model Fit: Observed vs Predicted",
        legend=:topleft,
        alpha=0.6,
        color=:blue,
        ms=4
    )

    # Add the predicted points (Stars)
    scatter!(sim_plot, trap_df.date, y_hat; 
             shape=:star, 
             color=:red, 
             ms=6,
             label="Model Prediction")

    # =========================================================
    # PLOT 2: RESIDUALS CHART
    # =========================================================

    resid_plot = scatter(
        trap_df.date, resid;
        label="Residual (Obs - Pred)",
        xlabel="Date",
        ylabel="Error",
        title="Residuals (Model Error)",
        legend=:topleft,
        color=:purple,
        markerstrokecolor=:white,
        alpha=0.8
    )

    hline!(resid_plot, [0.0]; label="Zero Error", lw=2, color=:black, linestyle=:dash)

    return (sim_plot=sim_plot, resid_plot=resid_plot, sim=sim)
end

"""
    summarize_mosquito_fit(trap_df, sim)

Calculates goodness-of-fit metrics (SSE, RMSE, MAE).
"""
function summarize_mosquito_fit(trap_df::DataFrame, sim)
    # Ensure column names match what is in trap_df
    if hasproperty(trap_df, :mfai_obvs)
        y_obs = Float64.(trap_df.mfai_obvs)
    elseif hasproperty(trap_df, :mfai)
        y_obs = Float64.(trap_df.mfai)
    else
        error("DataFrame missing observed data column (:mfai_obvs or :mfai)")
    end

    y_hat = sim.mfai_pred
    resid = y_obs .- y_hat

    # Calculate Metrics
    sse  = sum(resid .^ 2)           # Sum of Squared Errors
    rmse = sqrt(mean(resid .^ 2))    # Root Mean Squared Error
    mae  = mean(abs.(resid))         # Mean Absolute Error

    # --- PRINTING THE OUTPUT ---
    println("\n" * "="^30)
    println("   GOODNESS OF FIT METRICS   ")
    println("="^30)
    @printf("SSE (Sum Squared Error):  %.4f\n", sse)
    @printf("RMSE (Root Mean Sq Err):  %.4f\n", rmse)
    @printf("MAE (Mean Abs Error):     %.4f\n", mae)
    println("="^30 * "\n")

    return (SSE=sse, RMSE=rmse, MAE=mae)
end
"""
    mosquito_mfai_comp_table(trap_df, sim)

Returns a DataFrame with observed MFAI and predicted MFAI at each observed timestamp.
Also includes residuals (obs - pred).
"""
function mosquito_mfai_comp_table(trap_df::DataFrame, sim)
    # --- observed MFAI ---
    if hasproperty(trap_df, :mfai_obvs)
        mfai_observed = Float64.(trap_df.mfai_obvs)
    elseif hasproperty(trap_df, :mfai)
        mfai_observed = Float64.(trap_df.mfai)
    else
        error("trap_df missing observed MFAI column (:mfai_obvs or :mfai)")
    end

    # --- observed times (days since first trap date) ---
    start_date = minimum(trap_df.date)
    observed_times = Float64[Dates.value(d - start_date) for d in trap_df.date]

    # --- predicted MFAI ---
    # Prefer the precomputed vector from simulation output (what your plots use)
    mfai_predicted = if hasproperty(sim, :mfai_pred) && length(sim.mfai_pred) == length(observed_times)
        Float64.(sim.mfai_pred)
    else
        # Fallback: compute from the ODE solution at observed times
        # (works if `sim` is an ODESolution-like object)
        MosquitoModelDynamics.compute_mfai_theo(
            sim,
            observed_times,
            MosquitoModelDynamics.IX_T,
            Constants.N_TRAPS
        )
    end

    df = DataFrame(
        date = trap_df.date,
        t_sim_days = observed_times,
        mfai_observed = mfai_observed,
        mfai_predicted = mfai_predicted,
    )
    df.residual = df.mfai_observed .- df.mfai_predicted
    return df
end

end # module