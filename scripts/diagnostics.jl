using Printf
using DataInterpolations
#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using CSV
using DataFrames
using DataInterpolations
using Plots
using DifferentialEquations

# Import your package modules
using DengueODES
using DengueODES.Shared.Temperature  # This works now!
using DengueODES.Shared.TimeUtil
using DengueODES.MosquitoCapture
using DengueODES.MosquitoCapture.Report 



# 1. MOCK CONSTANTS (So the script runs standalone)
# In your real model, these come from Constants.jl
const N_HOUSEHOLDS = 107533.0 

# 2. DEFINE THE FUNCTIONS DIRECTLY
# (No module Entomology wrapper needed here)

@inline function heaviside(x)
    return x >= 0.0 ? 1.0 : 0.0
end

function get_carrying_capacity(t, C₀, b_cap, ϵ, t_start)
    dt = t - 6250 - ϵ
    u_t = heaviside(dt)
    final_capacity = C₀ + b_cap * dt * u_t
    return final_capacity * N_HOUSEHOLDS
end

function oviposition_rate(temp)
    b₀, b₁, b₂, b₃, b₄ = -5.3999, 1.800160, -2.12e-1, 1.02e-2, -1.51e-4
    δ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(0.0, δ)
end

function aquatic_transition(temp)
    if temp < 15.0 || temp > 35.0
        return 0.0
    else
        days = 21.0 - 0.4 * temp
        return 1.0 / max(days, 5.0)
    end
end

function aquatic_mortality(temp)
    b₀, b₁, b₂, b₃, b₄ = 2.130, -3.797e-1, 2.457e-2, -6.778e-4, 6.794e-6
    μ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(1e-6, μ)
end

function adult_mortality(temp)
    b₀, b₁, b₂, b₃, b₄ = 8.69e-1, -1.59e-1, 1.12e-2, -3.41e-4, 3.81e-6
    μₘ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(1e-6, μₘ)
end



# -------------------------------------------------------
# 1. Load Data (Trap + Weather)
# -------------------------------------------------------

# A. Process Weather & Create Interpolator
# This now calls your cleaned Temperature.jl

# B. Process Trap Data
# Assuming Data.get_trap_data() is the correct function name in your current setup
trap_df = DengueODES.MosquitoCapture.Data.get_trap_data()
sim_start_date = minimum(trap_df.date)
sim_end_date   = maximum(trap_df.date)

weather_df = Temperature.get_weather_data(force_process=true)

# C. CUT WEATHER DATA
# We keep 1 day before/after as a safety buffer for the solver
buffer = Day(1)
weather_subset = filter(row -> (sim_start_date - buffer) <= row.date <= (sim_end_date + buffer), weather_df)

# D. Create the Interpolator with the SUBSET
temp_interp = Temperature.get_temperature_interpolator(weather_subset)

# 3. DIAGNOSTIC EXECUTION
println("\n" * "="^80)
println("TIME SERIES DEBUGGING")
println("="^80)

# First: CHECK YOUR TEMPERATURE INTERPOLATOR
println("\n1. TEMPERATURE OVER TIME (every 90 days for 5 years)")
println("-"^60)
println("Day  | Year | Temp(°C) | Expected Pattern")
println("-"^60)

for t in 0:90:1825
    year = t/365
    T = temp_interp(t)  # YOUR INTERPOLATOR
    season = mod(t, 365) < 90 ? "Winter" : 
             mod(t, 365) < 180 ? "Spring" :
             mod(t, 365) < 270 ? "Summer" : "Fall"
    @printf("%4d | %.2f | %7.2f  | %s\n", t, year, T, season)
end

# Second: CHECK RATES OVER TIME
println("\n2. RATES OVER TIME (every 90 days)")
println("-"^80)
println("Day  | Temp | δ (ovip) | γₘ (emerg) | μₐ (aq.mort) | μₘ (ad.mort)")
println("-"^80)

for t in 0:90:1825
    T = temp_interp(t)
    δ = oviposition_rate(T)
    γm = aquatic_transition(T)
    μa = aquatic_mortality(T)
    μm = adult_mortality(T)
    @printf("%4d | %.1f | %8.4f | %10.6f | %12.6f | %12.6f\n", 
            t, T, δ, γm, μa, μm)
end

# Third: STATISTICAL CHECK
println("\n3. VARIABILITY CHECK")
println("-"^60)

temps = [temp_interp(t) for t in 0:10:1825]
δs = [oviposition_rate(temp_interp(t)) for t in 0:10:1825]
γms = [aquatic_transition(temp_interp(t)) for t in 0:10:1825]

println("Temperature: min=$(minimum(temps)), max=$(maximum(temps)), range=$(maximum(temps)-minimum(temps))")
println("Oviposition: min=$(minimum(δs)), max=$(maximum(δs)), range=$(maximum(δs)-minimum(δs))")
println("Emergence:   min=$(minimum(γms)), max=$(maximum(γms)), range=$(maximum(γms)-minimum(γms))")

# Fourth: LINEAR FIT CHECK
println("\n4. HOW LINEAR ARE YOUR RATES?")
println("(R² close to 1.0 = very linear, close to 0 = nonlinear)")
println("-"^60)

using Statistics
t_sample = collect(0:10:1825)
δ_sample = [oviposition_rate(temp_interp(t)) for t in t_sample]

# Simple linear regression
n = length(t_sample)
mean_t = mean(t_sample)
mean_δ = mean(δ_sample)
slope = sum((t_sample .- mean_t) .* (δ_sample .- mean_δ)) / sum((t_sample .- mean_t).^2)
intercept = mean_δ - slope * mean_t

# R²
ss_tot = sum((δ_sample .- mean_δ).^2)
ss_res = sum((δ_sample .- (slope .* t_sample .+ intercept)).^2)
r_squared = 1 - ss_res/ss_tot

println("Oviposition rate vs time: R² = $(round(r_squared, digits=4))")
println("  (If R² > 0.95, it's essentially linear)")
# println("="^80)
# println("MOSQUITO BIOLOGY DIAGNOSTIC SCRIPT")
# println("="^80)

# println("-"^80)
# @printf("%-10s | %-10s | %-10s | %-10s | %-10s\n", "Temp(°C)", "Oviposition", "Emergence", "Aq. Mort", "Ad. Mort")
# println("-"^80)

# for T in temps
#     δ = oviposition_rate(T)
#     γ = aquatic_transition(T)
#     μa = aquatic_mortality(T)
#     μm = adult_mortality(T)
#     @printf("%-10.1f | %-10.4f | %-10.4f | %-10.4f | %-10.4f\n", T, δ, γ, μa, μm)
# end

# println("-"^80)
# println("Run complete. Use these values to verify your model's temperature sensitivity.")