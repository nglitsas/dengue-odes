module Report

using DataFrames
using DataInterpolations
using Plots
using Statistics
using Dates
using CSV

# Import Shared utilities relatively
using ...Shared.TimeUtil

# Import Simulation relatively
using ..Simulate

# EXPORT ALL NECESSARY FUNCTIONS
export plot_mosquito_simulation, summarize_mosquito_fit, mosquito_mfai_comp_table

# ==========================================
# 1. HELPER: ALIGNMENT
# ==========================================
function align_simulation_to_observations(sim_df::DataFrame, trap_df::DataFrame)
    if nrow(sim_df) < 2
        # Return a safe fallback if simulation failed/empty
        return zeros(nrow(trap_df)) 
    end

    # Ensure simulation data is sorted and has no duplicates
    sort!(sim_df, :t)
    unique!(sim_df, :t)
    
    # Extract sorted arrays
    t_vals = Vector{Float64}(sim_df.t)
    mfai_vals = Vector{Float64}(sim_df.mfai_pred)
    
    # DEBUG: Print what we have
    println("\nDEBUG align_simulation_to_observations:")
    println("  nrow(sim_df) = ", nrow(sim_df))
    println("  First 5 t_vals: ", t_vals[1:min(5, length(t_vals))])
    println("  First 5 mfai_vals: ", mfai_vals[1:min(5, length(mfai_vals))])
    println("  Last 5 t_vals: ", t_vals[max(1,end-4):end])
    println("  Last 5 mfai_vals: ", mfai_vals[max(1,end-4):end])
    println("  issorted(t_vals) = ", issorted(t_vals))
    println("  issorted(mfai_vals) = ", issorted(mfai_vals))
    println("  any(isnan, t_vals) = ", any(isnan, t_vals))
    println("  any(isnan, mfai_vals) = ", any(isnan, mfai_vals))
    println("  length(t_vals) = ", length(t_vals))
    println("  length(mfai_vals) = ", length(mfai_vals))
    
    # Verify data integrity
    if !issorted(t_vals)
        error("Time values are not sorted even after sorting!")
    end
    
    if any(isnan, t_vals) || any(isnan, mfai_vals)
        error("NaN values found in interpolation data!")
    end
    
    if length(t_vals) != length(mfai_vals)
        error("Length mismatch: t_vals has $(length(t_vals)) elements, mfai_vals has $(length(mfai_vals))")
    end
    
    # Create interpolation: LinearInterpolation(x_data, y_data)
    # where x = time (independent variable), y = MFAI (dependent variable)
    println("  Creating LinearInterpolation with $(length(t_vals)) points...")
    
    interp = LinearInterpolation(mfai_vals, t_vals)  # DataInterpolations uses (u, t) order!
    
    println("  Interpolation created successfully!")
    
    # Evaluate at observation times
    pred_at_obs = Float64[]
    
    println("  Interpolating at $(length(trap_df.t)) observation times...")
    
    for (i, t) in enumerate(trap_df.t)
        # Handle extrapolation carefully
        if t < minimum(t_vals) || t > maximum(t_vals)
            # Use nearest neighbor for out-of-bounds
            if t < minimum(t_vals)
                push!(pred_at_obs, mfai_vals[1])
                if i <= 3
                    println("    t[$i]=$t is before sim range, using first value $(mfai_vals[1])")
                end
            else
                push!(pred_at_obs, mfai_vals[end])
                if i <= 3
                    println("    t[$i]=$t is after sim range, using last value $(mfai_vals[end])")
                end
            end
        else
            val = interp(t)
            push!(pred_at_obs, val)
            if i <= 3
                println("    t[$i]=$t → MFAI=$val")
            end
        end
    end
    
    println("  Interpolation complete! Generated $(length(pred_at_obs)) predictions\n")
    
    return pred_at_obs
end

# ==========================================
# 2. PLOTTING
# ==========================================
function plot_mosquito_simulation(trap_df::DataFrame, temp_interp; 
                                  fitted_params=nothing,
                                  solver=nothing, 
                                  saveat_daily=1.0,
                                  reltol=1e-6,
                                  abstol=1e-6)
    
    # Run the simulation
    results = Simulate.run_mosquito_simulation(
        trap_df, temp_interp; 
        fitted_params=fitted_params,
    )

    sim_df = results.sim
    rates_df = results.rates

    # Plot 1: Dynamics
    p1 = plot(sim_df.t, sim_df.mfai_pred, 
        label="Model (MFAI)", lw=2, color=:blue,
        xlabel="Time (t)", ylabel="Mosquito Abundance (Index)",
        title="Mosquito Population Dynamics")
    
    scatter!(p1, trap_df.t, trap_df.mfai_obvs, 
        label="Trap Data", color=:red, ms=3, alpha=0.6)

    # Plot 2: Residuals
    preds = align_simulation_to_observations(sim_df, trap_df)
    resids = preds .- trap_df.mfai_obvs
    
    p2 = scatter(trap_df.t, resids, 
        label="Residuals", color=:purple, ms=3,
        xlabel="Time (t)", ylabel="Model - Obs",
        title="Fit Residuals")
    hline!(p2, [0.0], color=:black, ls=:dash, label="")
    
    p_rates = plot_biological_rates(rates_df)

    return (
        sim = sim_df,
        sim_plot = p1,
        resid_plot = p2,
        rates_plot = p_rates
    )
    
end

# ==========================================
# 3. STATISTICS & TABLES
# ==========================================

function summarize_mosquito_fit(trap_df, sim_df)
    preds = align_simulation_to_observations(sim_df, trap_df)
    obs = trap_df.mfai_obvs
    
    # Handle edge case where simulation failed completely
    if all(preds .== 0.0) && mean(obs) > 0
        return (N=length(obs), MSE=Inf, SSE=Inf, RMSE=Inf, MAE=Inf, Correlation=0.0)
    end

    n = length(obs)
    resids = preds .- obs
    
    # --- ADDED SSE CALCULATION ---
    sse = sum(resids.^2)
    mse = mean(resids.^2)
    rmse = sqrt(mse)
    mae = mean(abs.(resids))
    
    r_val = cor(preds, obs)
    
    return (
        N = n,
        MSE = mse,
        SSE = sse,
        RMSE = rmse,
        MAE = mae,
        Correlation = r_val
    )
end

"""
    plot_biological_rates(df_rates)

Generates a multi-panel plot showing how climate drives mosquito biology.
"""
function plot_biological_rates(df_rates)
    dates = TimeUtil.t_to_date.(df_rates.t)

    # Plot 1a: Full range (log scale)
    p1a = plot(dates, df_rates.capacity, 
               title="Carrying Capacity (Log Scale)", 
               ylabel="Total Population", 
               label="Capacity", 
               color=:green, lw=2,
               yscale=:log10)
    
    # Plot 1b: Zoom on early times
    early_idx = df_rates.t .<= 7500
    p1b = plot(dates[early_idx], df_rates.capacity[early_idx],
               title="Carrying Capacity (Early, Linear)", 
               ylabel="Total Population", 
               label="Capacity", 
               color=:green, lw=2)

    # Plot 2a: Oviposition (SEPARATE)
    p2a = plot(dates, df_rates.oviposition_rate,
               title="Oviposition Rate", 
               ylabel="δ (eggs/female/day)", 
               label="Oviposition",
               color=:blue, lw=1.5)

    # Plot 2b: Emergence (SEPARATE) 
    p2b = plot(dates, df_rates.emergence_rate,
               title="Emergence Rate", 
               ylabel="γₘ (1/day)", 
               label="Emergence",
               color=:orange, lw=1.5)

    # Plot 3a: Aquatic Mortality (SEPARATE)
    p3a = plot(dates, df_rates.mu_a,
               title="Aquatic Mortality", 
               ylabel="μₐ (1/day)", 
               label="Aquatic",
               color=:red, lw=1.5)

    # Plot 3b: Adult Mortality (SEPARATE)
    p3b = plot(dates, df_rates.mu_m,
               title="Adult Mortality", 
               ylabel="μₘ (1/day)", 
               label="Adult",
               color=:purple, lw=1.5)

    return plot(p1a, p1b, p2a, p2b, p3a, p3b, 
                layout=(6,1), size=(800, 1800), xlabel="Date")
end

function mosquito_mfai_comp_table(trap_df::DataFrame, sim_df::DataFrame)
    # Align model predictions to the exact times of the trap observations
    preds = align_simulation_to_observations(sim_df, trap_df)
    
    # Create a clean comparison table
    comp_df = DataFrame(
        date = trap_df.date,
        t = trap_df.t,
        observed_mfai = trap_df.mfai_obvs,
        predicted_mfai = preds,
        residual = preds .- trap_df.mfai_obvs
    )
    
    return comp_df
end

end # module