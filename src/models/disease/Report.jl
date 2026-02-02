module Report

using DataFrames
using Plots
using Dates
using Statistics
using Printf

# Import siblings
using ..Simulate
using ..Dengue # Access to ModelParams if needed

export plot_dengue_simulation, summarize_dengue_fit

"""
    plot_dengue_simulation(case_df, temp_interp; reporting_rate=1.0, ...)

Runs the simulation and plots Observed Cases vs Predicted Incidence.
"""
function plot_dengue_simulation(case_df::DataFrame,
                                temp_interp;
                                reporting_rate::Float64 = 1.0, # Scale down model to match data
                                fitted_params = nothing,
                                solver = Rosenbrock23())

    # 1. Run Simulation
    sim_results = Simulate.run_dengue_simulation(
        case_df,
        temp_interp;
        fitted_params = fitted_params,
        solver = solver
    )

    # 2. Extract Data
    # Apply reporting rate (Model predicts All Infections -> Data is Reported Cases)
    pred_cases = sim_results.incidence .* reporting_rate
    
    # 3. Create Plot
    p1 = plot(
        case_df.dt_sin_pri, 
        case_df.notified, 
        seriestype = :scatter, 
        label = "Observed Cases (Notified)",
        color = :black,
        alpha = 0.6,
        markersize = 3,
        title = "Dengue Model Fit",
        ylabel = "Daily Cases",
        xlabel = "Date"
    )

    plot!(
        p1,
        sim_results.dates,
        pred_cases,
        linewidth = 2,
        color = :red,
        label = "Predicted (Rate=$(reporting_rate))"
    )

    # 4. Create Residual Plot
    # We need to interpolate model predictions to the exact days data exists
    # (Though usually they are both daily, so we can just match dates)
    
    # Simple matching (assuming daily data and daily saveat)
    # Be careful: Simulation might have 1 extra day or slight mismatch
    
    # Filter simulation to match data dates
    common_dates = intersect(case_df.dt_sin_pri, sim_results.dates)
    
    # Filter data to common dates
    obs_subset = filter(row -> row.dt_sin_pri in common_dates, case_df)
    
    # Find indices in simulation that match these dates
    sim_indices = [findfirst(==(d), sim_results.dates) for d in obs_subset.dt_sin_pri]
    pred_subset = pred_cases[sim_indices]
    
    residuals = obs_subset.notified .- pred_subset

    p2 = plot(
        obs_subset.dt_sin_pri,
        residuals,
        seriestype = :bar,
        color = :blue,
        alpha = 0.5,
        title = "Residuals (Observed - Predicted)",
        ylabel = "Error",
        legend = false
    )

    # Combine
    final_plot = plot(p1, p2, layout=(2,1), size=(800, 800))

    return (
        plot = final_plot,
        sim = sim_results,
        residuals = residuals
    )
end

function summarize_dengue_fit(case_df, sim_results, reporting_rate=1.0)
    # Calculate simple error metrics
    # Logic similar to residuals block above...
    common_dates = intersect(case_df.dt_sin_pri, sim_results.dates)
    obs_subset = filter(row -> row.dt_sin_pri in common_dates, case_df)
    sim_indices = [findfirst(==(d), sim_results.dates) for d in obs_subset.dt_sin_pri]
    
    obs = Float64.(obs_subset.notified)
    pred = sim_results.incidence[sim_indices] .* reporting_rate

    sse  = sum((obs .- pred).^2)
    mse  = sse / length(obs)
    rmse = sqrt(mse)
    mae  = mean(abs.(obs .- pred))

    return (
        N    = length(obs),
        SSE  = sse,           
        RMSE = rmse,
        MAE  = mae
    )
end

end # module