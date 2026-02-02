#!/usr/bin/env julia
"""
COMPREHENSIVE DEBUG SCRIPT

This will help us figure out why the model predicts zero mosquitoes.
Add this to your scripts folder and run it.
"""

using Pkg
Pkg.activate(".")

using DifferentialEquations
using DataFrames
using CSV
using Plots
using Printf

# Load your modules
using DengueODES
using DengueODES.Shared.TimeUtil
using DengueODES.Shared.Temperature
using DengueODES.MosquitoCapture.Data
using DengueODES.MosquitoCapture.Constants
using DengueODES.MosquitoCapture.Entomology
using DengueODES.MosquitoCapture.MosquitoModelDynamics

println("="^70)
println("COMPREHENSIVE MODEL DEBUG")
println("="^70)

# Load data
trap_df = Data.get_trap_data()
weather_df = Temperature.get_weather_data(force_process=false)
temp_interp = Temperature.get_temperature_interpolator(weather_df)

# Use fitted parameters
C₀ = 300000.0
b_cap = 0.349
ε = 974.0

# Get first trap time
t_start = TimeUtil.date_to_t(trap_df.date[1])
println("\nData Information:")
println("  First trap date: $(trap_df.date[1])")
println("  t_start = $t_start")
println("  C₀ = $C₀")
println("  b_cap = $b_cap")
println("  ε = $ε")

# Create model parameters
p = MosquitoModelParams(C₀, b_cap, ε, t_start, temp_interp)

# Get initial conditions
u0 = default_u0(p)
println("\nInitial Conditions:")
println("  A(0) = $(u0[1])")
println("  M(0) = $(u0[2])")
println("  T(0) = $(u0[3])")
println("  Ratio M/A = $(u0[2]/u0[1])")

# Test carrying capacity
println("\nCarrying Capacity Test:")
for t_offset in [0.0, 100.0, 500.0, 900.0, 1000.0, 1500.0]
    t = t_start + t_offset
    C = Entomology.get_carrying_capacity(t, C₀, b_cap, ε, t_start)
    println("  C(t=$(round(Int, t))) = $(round(Int, C))")
end

# Test temperature and rates
println("\nTemperature and Rates at t_start:")
temp = temp_interp(t_start)
println("  Temperature = $(round(temp, digits=2))°C")

δₜ, γₘₜ, μₐₜ, μₘₜ = Entomology.get_rates(t_start, temp_interp)
println("  Oviposition rate (δ) = $(round(δₜ, digits=4))")
println("  Aquatic transition (γₘ) = $(round(γₘₜ, digits=4))")
println("  Aquatic mortality (μₐ) = $(round(μₐₜ, digits=4))")
println("  Adult mortality (μₘ) = $(round(μₘₜ, digits=4))")

# Test trapping rate
trap_rate = Constants.ALPHA * (Constants.N_TRAPS / Constants.N_HOUSEHOLDS)
println("\nTrapping:")
println("  α = $(Constants.ALPHA)")
println("  Ntraps = $(Constants.N_TRAPS)")
println("  Nhouseholds = $(Constants.N_HOUSEHOLDS)")
println("  Trap rate = $(round(trap_rate, digits=6))")
println("  Expected trapped per day = $(round(trap_rate * u0[2], digits=1))")

# Run a short simulation
println("\n" * "="^70)
println("RUNNING SHORT SIMULATION (10 days)")
println("="^70)

tspan = (t_start, t_start + 10.0)
prob = ODEProblem(CaptureModel_Fast, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=1.0)

println("\nSimulation Results:")
println("  Success: $(sol.retcode == ReturnCode.Success)")
println("\nState variables over time:")
println("  Time    A          M          T")
println("  " * "-"^50)

for i in 1:min(11, length(sol.t))
    t = sol.t[i]
    A, M, T = sol.u[i]
    @printf("  %.1f   %.0f   %.0f   %.2f\n", t, A, M, T)
end

# Check if mosquitoes are actually being trapped
println("\nTrapped Mosquitoes Analysis:")
T_initial = sol.u[1][3]
T_final = sol.u[end][3]
T_delta = T_final - T_initial

println("  T(0) = $(round(T_initial, digits=2))")
println("  T(10) = $(round(T_final, digits=2))")
println("  Delta = $(round(T_delta, digits=2))")
println("  MFAI per 10 days = $(round(T_delta / Constants.N_TRAPS, digits=4))")

if T_delta < 1.0
    println("\n⚠️  CRITICAL ISSUE: Almost no mosquitoes trapped!")
    println("     Possible causes:")
    println("     1. M(t) dropping to zero quickly")
    println("     2. Trap rate calculation wrong")
    println("     3. dT/dt equation not working")
    println("     4. Numerical issues")
end

# Check if M is crashing
M_initial = sol.u[1][2]
M_final = sol.u[end][2]
M_ratio = M_final / M_initial

println("\nAdult Mosquito Population Analysis:")
println("  M(0) = $(round(M_initial, digits=0))")
println("  M(10) = $(round(M_final, digits=0))")
println("  Ratio = $(round(M_ratio, digits=3))")

if M_ratio < 0.5
    println("\n⚠️  WARNING: M dropped by >50% in 10 days!")
    println("     This suggests mortality is too high or emergence too low")
end

# Now test at actual collection times
println("\n" * "="^70)
println("TESTING AT COLLECTION TIMES")
println("="^70)

collection_times = Float64.(TimeUtil.date_to_t.(trap_df.date[1:5]))
tspan_full = (collection_times[1], collection_times[end])

prob_full = ODEProblem(CaptureModel_Fast, u0, tspan_full, p)
sol_full = solve(prob_full, Tsit5(), saveat=collection_times)

println("\nCollection Results:")
println("  Collection    t        T_cumul    Delta     MFAI_calc    MFAI_obs")
println("  " * "-"^70)

for i in 1:5
    t = collection_times[i]
    T_curr = sol_full(t)[3]
    T_prev = i == 1 ? 0.0 : sol_full(collection_times[i-1])[3]
    delta = T_curr - T_prev
    mfai_calc = delta / Constants.N_TRAPS
    mfai_obs = trap_df.mfai_obvs[i]
    
    @printf("  %d         %.1f    %.2f       %.2f      %.4f       %.4f\n", 
            i, t, T_curr, delta, mfai_calc, mfai_obs)
end

println("\n" * "="^70)
println("DIAGNOSIS")
println("="^70)

# Final diagnosis
if T_delta < 1.0
    println("\n❌ MODEL IS NOT TRAPPING MOSQUITOES")
    println("\nLikely causes (in order of probability):")
    println("1. M(t) crashes to near zero immediately")
    println("   → Check temperature-dependent mortality rates")
    println("   → Check if A(0)/C(0) ratio allows oviposition")
    println("2. Trap rate is being calculated incorrectly")
    println("   → Verify dT/dt = trap_rate * M in ODE")
    println("3. Initial conditions are wrong")
    println("   → Verify u0 is being used correctly")
elseif T_delta > 100
    println("\n✅ MODEL IS TRAPPING MOSQUITOES!")
    println("\nBut MFAI is still wrong. Check:")
    println("1. Are you dividing by the right number of traps?")
    println("2. Are collection times aligned correctly?")
    println("3. Is there a units mismatch somewhere?")
end

println("\n" * "="^70)
println("END DEBUG")
println("="^70)
