#!/usr/bin/env julia
"""
DIAGNOSTIC SCRIPT: Check Model Output Scale

This script will help you understand why your model might not be fitting the data.
Run this AFTER you've made the corrections to see what's happening.
"""

using CSV
using DataFrames
using Plots
using Printf

println("="^70)
println("MOSQUITO MODEL OUTPUT DIAGNOSTICS")
println("="^70)

# Load your observed data
trap_df = CSV.read("data/mosq_trapped_model.csv", DataFrame)

println("\n1. OBSERVED DATA ANALYSIS:")
println("-"^70)
println("Date range: $(minimum(trap_df.date)) to $(maximum(trap_df.date))")
println("Number of collections: $(nrow(trap_df))")
println()
println("MFAI statistics:")
println("  Min:    $(round(minimum(trap_df.mfai_obvs), digits=4))")
println("  Max:    $(round(maximum(trap_df.mfai_obvs), digits=4))")
println("  Mean:   $(round(mean(trap_df.mfai_obvs), digits=4))")
println("  Median: $(round(median(trap_df.mfai_obvs), digits=4))")
println()

println("Trapped mosquitoes statistics:")
println("  Min:    $(minimum(trap_df.trapped))")
println("  Max:    $(maximum(trap_df.trapped))")
println("  Mean:   $(round(mean(trap_df.trapped), digits=1))")
println("  Median: $(round(median(trap_df.trapped), digits=1))")
println()

println("Number of traps statistics:")
println("  Min:    $(minimum(trap_df.n_traps))")
println("  Max:    $(maximum(trap_df.n_traps))")
println("  Mean:   $(round(mean(trap_df.n_traps), digits=1))")
println()

# Verify MFAI calculation
println("2. VERIFY MFAI CALCULATION:")
println("-"^70)
for i in 1:min(5, nrow(trap_df))
    trapped = trap_df.trapped[i]
    ntraps = trap_df.n_traps[i]
    mfai_obs = trap_df.mfai_obvs[i]
    mfai_calc = trapped / ntraps
    
    match = abs(mfai_calc - mfai_obs) < 1e-6 ? "✓" : "✗"
    
    println("Row $i: Trapped=$trapped, Ntraps=$ntraps")
    println("  Observed MFAI: $(round(mfai_obs, digits=4))")
    println("  Calculated:    $(round(mfai_calc, digits=4)) $match")
    println()
end

println("\n3. EXPECTED MODEL BEHAVIOR:")
println("-"^70)

# Paper parameters
N = 256088.0
N_TRAPS = 2412.0
N_HOUSEHOLDS = 102751.0
ALPHA = 0.02
K = 0.5

# Paper's fitted values (interpreted as density)
C0_density = 1.27  # per household
bcap = 0.3165
epsilon = 909.0

# Total capacity
C0_total = C0_density * N_HOUSEHOLDS

println("Paper's parameters:")
println("  C₀ (density):  $C0_density mosquitoes/household")
println("  C₀ (total):    $(round(Int, C0_total)) mosquitoes")
println("  b_cap:         $bcap")
println("  ε:             $epsilon days")
println()

# Initial conditions
A0 = 0.85 * C0_total
M0 = 0.7 * N
println("Initial conditions:")
println("  A(0) = $(round(Int, A0)) aquatic mosquitoes")
println("  M(0) = $(round(Int, M0)) adult mosquitoes")
println()

# Trap rate
trap_rate = ALPHA * (N_TRAPS / N_HOUSEHOLDS)
println("Trapping dynamics:")
println("  α = $ALPHA")
println("  Ntraps/Nhouseholds = $(round(N_TRAPS/N_HOUSEHOLDS, digits=4))")
println("  Trap rate = α × (Ntraps/Nhouseholds) = $(round(trap_rate, sigdigits=4))")
println()

# Expected catches
println("Expected catches over 60 days (typical collection interval):")
println()

# Assume M stays roughly constant at M0
days = 60.0
trapped_per_60d = trap_rate * M0 * days
mfai_per_60d = trapped_per_60d / N_TRAPS

println("  If M(t) ≈ M(0) = $(round(Int, M0)):")
println("    Trapped in 60 days: $(round(Int, trapped_per_60d))")
println("    MFAI: $(round(mfai_per_60d, digits=2))")
println()

# Compare with observed
median_trapped = median(trap_df.trapped)
median_mfai = median(trap_df.mfai_obvs)

println("  Observed median:")
println("    Trapped in collection: $(round(Int, median_trapped))")
println("    MFAI: $(round(median_mfai, digits=2))")
println()

ratio = mfai_per_60d / median_mfai
println("  Ratio: Model/Observed = $(round(ratio, digits=1))×")
println()

if ratio > 5
    println("  ⚠️  MODEL IS TRAPPING $(round(ratio, digits=1))× TOO MANY MOSQUITOES!")
    println()
    println("  Possible causes:")
    println("    1. α = 0.02 might be too high")
    println("    2. M(t) might be staying too high")
    println("    3. Trap rate calculation might be wrong")
    println("    4. Initial M(0) might be too high")
elseif ratio < 0.2
    println("  ⚠️  MODEL IS TRAPPING $(round(1/ratio, digits=1))× TOO FEW MOSQUITOES!")
else
    println("  ✓ Model is in the right ballpark!")
end

println("\n4. PARAMETER SEARCH SPACE:")
println("-"^70)

println("Your current search space:")
println("  C₀: [0.1, 5.0] (density per household)")
println("  b_cap: [0.0, 1.2]")
println("  ε: [0, 365] days")
println()

println("Recommended search space (to match paper):")
println("  C₀: [0.5, 3.0] (density per household)")
println("  b_cap: [0.0, 1.2]")
println("  ε: [500, 1500] days")
println()

println("Paper's fitted values in your units:")
println("  C₀ ≈ 1.27 (density)")
println("  b_cap ≈ 0.3165")
println("  ε ≈ 909 days")
println()

println("\n5. SCALING ANALYSIS:")
println("-"^70)

# Check what happens with different M values
println("If we adjust M(t) or α to match observed MFAI:")
println()

target_mfai = median_mfai
target_trapped = target_mfai * N_TRAPS

println("Target: Trap $(round(Int, target_trapped)) mosquitoes in 60 days")
println("        to get MFAI = $(round(target_mfai, digits=3))")
println()

# What M is needed?
M_needed = target_trapped / (trap_rate * days)
println("Required M(t) with α=0.02:")
println("  M = $(round(Int, M_needed))")
println("  This is $(round(M_needed/M0, digits=2))× smaller than M(0)")
println()

# What alpha is needed?
alpha_needed = target_trapped / (M0 * days * (N_TRAPS/N_HOUSEHOLDS))
println("Required α with M=M(0):")
println("  α = $(round(alpha_needed, sigdigits=3))")
println("  This is $(round(alpha_needed/ALPHA, digits=2))× smaller than 0.02")
println()

println("="^70)
println("RECOMMENDATIONS:")
println("="^70)
println()
println("1. Check your model output: What is M(t) doing over time?")
println("   - Plot M(t) from your simulation")
println("   - It should vary between ~50k and ~300k")
println()
println("2. Verify trap rate calculation in your ODE:")
println("   - Should be: α × (Ntraps/Nhouseholds) × M(t)")
println("   - Check that constants match the paper")
println()
println("3. Consider whether α needs to be fitted:")
println("   - Paper assumes α = 0.02")
println("   - But this might need calibration")
println()
println("4. Update parameter search bounds:")
println("   - Use C₀ ∈ [0.5, 3.0]")
println("   - Use ε ∈ [500, 1500]")
println("="^70)

# Create a visualization
println("\nGenerating plot...")

p = plot(trap_df.date, trap_df.mfai_obvs,
         marker=:circle,
         markersize=6,
         label="Observed MFAI",
         xlabel="Date",
         ylabel="MFAI (mosquitoes per trap)",
         title="Observed MFAI Over Time",
         legend=:topleft,
         linewidth=2)

# Add horizontal lines for reference
hline!([median_mfai], label="Median", linestyle=:dash, color=:red)

savefig(p, "mosquito_diagnostics_plot.png")
println("Plot saved as 'mosquito_diagnostics_plot.png'")

println("\nDiagnostics complete!")
