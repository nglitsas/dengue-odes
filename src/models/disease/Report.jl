module Report

using DataFrames
using Plots
using Dates
using Statistics
using Printf

# Import siblings
using ...Shared.TimeUtil
using ..Simulate
using ..Dengue # Access to ModelParams if needed
using ..DengueModel


export plot_dengue_simulation, summarize_dengue_fit, plot_temperature_overlay,
 plot_prediction_scatter, plot_monthly_comparison,
       plot_human_compartments, plot_mosquito_compartments


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
# ============================================
# NEW PLOTTING FUNCTIONS
# ============================================

"""
Plot dengue cases with temperature overlay on secondary axis
"""
function plot_temperature_overlay(sim_results, cases_df, temp_interp, reporting_rate; 
                                  date_col=:dt_sin_pri)
    # Get temperature for all simulation dates
    temp_vals = [temp_interp(TimeUtil.date_to_t(d)) for d in sim_results.dates]
    
    p = plot(
        sim_results.dates, 
        sim_results.incidence .* reporting_rate,
        label = "Predicted Cases",
        linewidth = 2,
        color = :blue,
        ylabel = "Daily Cases",
        legend = :topleft,
        size = (1000, 500)
    )
    
    scatter!(p, 
        cases_df[!, date_col], 
        cases_df.notified,
        label = "Observed Cases",
        color = :red,
        markersize = 4,
        alpha = 0.6
    )
    
    p_temp = twinx(p)
    plot!(p_temp,
        sim_results.dates,
        temp_vals,
        label = "Temperature",
        linewidth = 2,
        color = :orange,
        ylabel = "Temperature (°C)",
        legend = :topright,
        linestyle = :dash
    )
    
    title!(p, "Dengue Cases vs Temperature")
    
    return p
end

"""
Create predicted vs observed scatter plot with R²
"""
function plot_prediction_scatter(sim_results, cases_df, reporting_rate; 
                                 date_col=:dt_sin_pri)
    # Find common dates
    common_dates = intersect(cases_df[!, date_col], sim_results.dates)
    
    if length(common_dates) == 0
        @warn "No overlapping dates for scatter plot"
        return nothing
    end
    
    obs_subset = filter(row -> row[date_col] in common_dates, cases_df)
    sim_indices = [findfirst(==(d), sim_results.dates) for d in obs_subset[!, date_col]]
    
    pred_aligned = sim_results.incidence[sim_indices] .* reporting_rate
    obs_aligned = obs_subset.notified
    
    # Calculate R²
    ss_res = sum((obs_aligned .- pred_aligned).^2)
    ss_tot = sum((obs_aligned .- mean(obs_aligned)).^2)
    r_squared = 1 - ss_res/ss_tot
    
    p = scatter(
        obs_aligned,
        pred_aligned,
        xlabel = "Observed Cases",
        ylabel = "Predicted Cases",
        title = "Predicted vs Observed\nR² = $(round(r_squared, digits=3))",
        label = "",
        markersize = 6,
        alpha = 0.6,
        size = (600, 600)
    )
    
    max_val = max(maximum(obs_aligned), maximum(pred_aligned))
    plot!(p, [0, max_val], [0, max_val], 
          label = "Perfect Prediction", 
          linewidth = 2, 
          linestyle = :dash,
          color = :red)
    
    return p
end

"""
Plot monthly aggregated cases comparison
"""
function plot_monthly_comparison(sim_results, cases_df, reporting_rate; 
                                 date_col=:dt_sin_pri)
    # Create comparison dataframe
    comparison_df = DataFrame(
        date = sim_results.dates,
        predicted_cases = sim_results.incidence .* reporting_rate
    )
    
    # Add observed data
    observed_dict = Dict(cases_df[!, date_col] .=> cases_df.notified)
    comparison_df[!, :observed_cases] = [get(observed_dict, d, missing) for d in comparison_df.date]
    
    # Aggregate by month
    comparison_df[!, :month] = Dates.yearmonth.(comparison_df.date)
    
    monthly_agg = combine(
        groupby(comparison_df, :month),
        :predicted_cases => sum => :predicted_total,
        :observed_cases => (x -> sum(skipmissing(x))) => :observed_total
    )
    
    # Convert to dates for plotting
    monthly_agg[!, :date] = [Date(ym[1], ym[2], 1) for ym in monthly_agg.month]
    
    p = plot(
        monthly_agg.date,
        monthly_agg.predicted_total,
        label = "Predicted",
        linewidth = 3,
        marker = :circle,
        markersize = 6,
        xlabel = "Month",
        ylabel = "Total Cases",
        title = "Monthly Case Comparison",
        size = (1000, 500),
        legend = :topleft
    )
    
    plot!(p,
        monthly_agg.date,
        monthly_agg.observed_total,
        label = "Observed",
        linewidth = 3,
        marker = :square,
        markersize = 6
    )
    
    return p
end

"""
Plot human compartment dynamics (SEIR)
"""
function plot_human_compartments(sim_results)
    if !all(hasfield(typeof(sim_results), f) for f in [:Hs, :He, :Hi, :Hr])
        @warn "Simulation results missing compartment data"
        return nothing
    end
    
    p = plot(
        sim_results.dates,
        [sim_results.Hs sim_results.He sim_results.Hi sim_results.Hr],
        label = ["Susceptible" "Exposed" "Infected" "Recovered"],
        linewidth = 2,
        xlabel = "Date",
        ylabel = "Number of Individuals",
        title = "Human Compartment Dynamics",
        legend = :right,
        size = (1000, 500)
    )
    
    return p
end

"""
Plot mosquito population dynamics
"""
function plot_mosquito_compartments(sim_results)
    if !all(hasfield(typeof(sim_results), f) for f in [:Ma, :Ms, :Me, :Mi])
        @warn "Simulation results missing mosquito compartment data"
        return nothing
    end
    
    p = plot(
        sim_results.dates,
        sim_results.Ma,
        label = "Aquatic (Ma)",
        linewidth = 2,
        xlabel = "Date",
        ylabel = "Mosquito Count",
        title = "Mosquito Population Dynamics",
        legend = :topright,
        size = (1000, 500)
    )
    
    plot!(p, sim_results.dates, sim_results.Ms, label = "Susceptible (Ms)", linewidth = 2)
    plot!(p, sim_results.dates, sim_results.Me, label = "Exposed (Me)", linewidth = 2)
    plot!(p, sim_results.dates, sim_results.Mi, label = "Infected (Mi)", linewidth = 2)
    
    return p
end


end # module