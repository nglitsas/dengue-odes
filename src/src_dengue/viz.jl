# src/viz.jl
# This file sets up the functions to visualize the model output
using Plots
import ..DengueTransmission: simulate, ModelParams, IX_Hₑ

function plot_fit(df, p::ModelParams; title_str="Dengue fit")
    t_obs = Float64.(df.t)
    y_obs = Float64.(df.cases)
    tspan = (minimum(t_obs), maximum(t_obs))

    sol = simulate(p, tspan; saveat=t_obs)
    H_e = Array(sol)[IX_Hₑ, :]
    yhat = p.θₕ .* H_e

    plt = plot(t_obs, y_obs; seriestype=:scatter, label="Observed cases",
               title=title_str, xlabel="t", ylabel="cases/incidence")
    plot!(plt, t_obs, yhat; label="Model θₕ·Hₑ (incidence proxy)")
    return plt
end