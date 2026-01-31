# scripts/run_pipeline.jl
# This file imports the data, calls the required functions, plots the fit, and outputs important statistics
using Pkg
Pkg.activate(@__DIR__ * "/..")

using DengueTransmissionProject
using Plots

df = load_cases_csv("data/raw_cases.csv")

base = ModelParams(
    N=1_000_000.0, Nⱼ=100.0, H₀=200_000.0, j=0.1,
    μₕ=1/(70*365), θₕ=1/5, γₕ=1/7,
    k=0.5, cₐ=0.0, cₘ=0.0,
    δₜ=10.0, μₘₜ=1/10, μₐₜ=1/7, γₘₜ=1/7, θₘₜ=1/10,
    bₜ=0.3, βₘₜ=0.3, βₕₜ=0.3,
    aᵦ=1.0, αₘ=1.0, αₕ=1.0,
    C=200_000.0, C₀=200_000.0, bₖ=1.0, ϵ=0.1,
    ϕ=0.99
)

bounds = Dict(
  :βₕₜ => (0.01, 2.0),
  :βₘₜ => (0.01, 2.0),
  :bₜ  => (0.05, 1.5),
  :θₘₜ => (0.05, 1.0)
)

best_p, best_loss = fit_lhs(df, base; bounds=bounds, nsamples=2000, seed=42)
println("best_loss = ", best_loss)
println("best params: βₕₜ=$(best_p.βₕₜ), βₘₜ=$(best_p.βₘₜ), bₜ=$(best_p.bₜ), θₘₜ=$(best_p.θₘₜ)")

plt = plot_fit(df, best_p; title_str="Dengue LHS Best Fit")
display(plt)