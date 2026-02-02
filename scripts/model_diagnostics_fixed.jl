#!/usr/bin/env julia
"""
DIAGNOSTIC SCRIPT: Check Model Output Scale

This script analyzes your mosquito trap data and compares it with expected model behavior.
Adapted to work with your project structure.
"""

using Pkg
Pkg.activate(".")  # Activate the current project

using CSV
using DataFrames
using Dates
using Statistics
using Plots
using Printf

println("="^70)
println("MOSQUITO MODEL OUTPUT DIAGNOSTICS")
println("="^70)

# Try to load processed data first, fall back to raw data
trap_df = nothing
processed_path = "data/processed/mosq_trapped_model.csv"
raw_path = "data/raw/mosq_aaeg_trap-2017_2022.csv"

if isfile(processed_path)
    println("\nLoading processed trap data from: $processed_path")
    trap_df = CSV.read(processed_path, DataFrame)
elseif isfile(raw_path)
    println("\nProcessed data not found. Loading raw data from: $raw_path")
    raw_df = CSV.read(raw_path, DataFrame)
    
    # Aggregate to monthly data (odd months only, similar to paper)
    println("  Raw data has $(nrow(raw_df)) daily records")
    println("  Need to aggregate to bimonthly collections...")
    println()
    println("  Note: Your raw data needs to be processed first.")
    println("  Please run your data processing pipeline to create bimonthly aggregates.")
    println()
    error("Cannot proceed without bimonthly aggregated data")
else
    error("Could not find trap data at:\n  $processed_path\n  $raw_path")
end

println("\n1. OBSERVED DATA ANALYSIS:")
println("-"^70)
println("Date range: $(minimum(trap_df.date)) to $(maximum(trap_df.date))")
println("Number of collections: $(nrow(trap_df))")
println()

# Find MFAI column (could be mfai_obvs or mfai)
mfai_col = nothing
if hasproperty(trap_df, :mfai_obvs)
    mfai_col = :mfai_obvs
elseif hasproperty(trap_df, :mfai)
    mfai_col = :mfai
else
    error("Could not find MFAI column (expected :mfai_obvs or :mfai)")
end

mfai_vals = trap_df[!, mfai_col]

println("MFAI statistics:")
println("  Min:    $(round(minimum(mfai_vals), digits=4))")
println("  Max:    $(round(maximum(mfai_vals), digits=4))")
println("  Mean:   $(round(mean(mfai_vals), digits=4))")
println("  Median: $(round(median(mfai_vals), digits=4))")
println()

if hasproperty(trap_df, :trapped)
    println("Trapped mosquitoes statistics:")
    println("  Min:    $(minimum(trap_df.trapped))")
    println("  Max:    $(maximum(trap_df.trapped))")
    println("  Mean:   $(round(mean(trap_df.trapped), digits=1))")
    println("  Median: $(round(median(trap_df.trapped), digits=1))")
    println()
end

if hasproperty(trap_df, :n_traps)
    println("Number of traps statistics:")
    println("  Min:    $(minimum(trap_df.n_traps))")
    println("  Max:    $(maximum(trap_df.n_traps))")
    println("  Mean:   $(round(mean(trap_df.n_traps), digits=1))")
    println()
end

# Verify MFAI calculation if we have trapped and n_traps
if hasproperty(trap_df, :trapped) && hasproperty(trap_df, :n_traps)
    println("2. VERIFY MFAI CALCULATION:")
    println("-"^70)
    for i in 1:min(5, nrow(trap_df))
        trapped = trap_df.trapped[i]
        ntraps = trap_df.n_traps[i]
        mfai_obs = mfai_vals[i]
        mfai_calc = trapped / ntraps
        
        match = abs(mfai_calc - mfai_obs) < 1e-6 ? "✓" : "✗"
        
        @printf("Row %d: Trapped=%d, Ntraps=%d\n", i, trapped, ntraps)
        @printf("  Observed MFAI: %.4f\n", mfai_obs)
        @printf("  Calculated:    %.4f %s\n\n", mfai_calc, match)
    end
end

println("\n3. EXPECTED MODEL BEHAVIOR:")
println("-"^70)

# Paper parameters
N = 256088.0
N_TRAPS = 2412.0
N_HOUSEHOLDS = 102751.0
ALPHA = 0.02
K = 0.5

# Paper's fitted values - treating C₀ as TOTAL capacity (not per-household)
C0_total = 1.27e5  # Total mosquitoes
bcap = 0.3165
epsilon = 909.0

println("Paper's parameters:")
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
println("  Trap rate per mosquito per day = $(round(trap_rate, sigdigits=4))")
println()

# Expected catches over 60 days
println("Expected catches over 60 days (typical collection interval):")
println()

days = 60.0
# Assuming M stays roughly constant at M0
trapped_per_60d_constant = trap_rate * M0 * days
mfai_per_60d_constant = trapped_per_60d_constant / N_TRAPS

println("  Scenario 1: If M(t) stays constant at M(0):")
println("    Daily trap rate: $(round(trap_rate * M0, digits=1)) mosquitoes/day")
println("    Trapped in 60 days: $(round(Int, trapped_per_60d_constant))")
println("    MFAI: $(round(mfai_per_60d_constant, digits=2))")
println()

# More realistic: M drops due to trapping
# dM/dt includes -trap_rate*M term, so M decreases exponentially
# M(t) = M0 * exp(-trap_rate * t)
# Trapped = integral from 0 to t of trap_rate * M0 * exp(-trap_rate*s) ds
#         = M0 * (1 - exp(-trap_rate * t))
trapped_per_60d_exp = M0 * (1 - exp(-trap_rate * days))
mfai_per_60d_exp = trapped_per_60d_exp / N_TRAPS

println("  Scenario 2: If M(t) decreases exponentially (more realistic):")
println("    M(60) = $(round(Int, M0 * exp(-trap_rate * days)))")
println("    Trapped in 60 days: $(round(Int, trapped_per_60d_exp))")
println("    MFAI: $(round(mfai_per_60d_exp, digits=2))")
println()

# Compare with observed
median_mfai = median(mfai_vals)

if hasproperty(trap_df, :trapped)
    median_trapped = median(trap_df.trapped)
    println("  Observed median:")
    println("    Trapped in collection: $(round(Int, median_trapped))")
    println("    MFAI: $(round(median_mfai, digits=3))")
    println()
    
    ratio_constant = mfai_per_60d_constant / median_mfai
    ratio_exp = mfai_per_60d_exp / median_mfai
    
    println("  Model/Observed ratios:")
    println("    Scenario 1 (constant M): $(round(ratio_constant, digits=1))×")
    println("    Scenario 2 (exponential): $(round(ratio_exp, digits=1))×")
    println()
    
    if ratio_exp > 5
        println("  ⚠️  MODEL IS TRAPPING $(round(ratio_exp, digits=1))× TOO MANY MOSQUITOES!")
        println()
        println("  Possible causes:")
        println("    1. α = 0.02 is too high for this trap system")
        println("    2. M(0) = 0.7N might be too high")
        println("    3. Temperature-dependent mortality might be too low")
        println("    4. Mosquito emergence rate might be too high")
    elseif ratio_exp < 0.2
        println("  ⚠️  MODEL IS TRAPPING $(round(1/ratio_exp, digits=1))× TOO FEW MOSQUITOES!")
    else
        println("  ✓ Model is in the right ballpark!")
    end
end

println("\n4. PARAMETER SEARCH SPACE:")
println("-"^70)

println("Your current search space in run_mosquito_fit.jl:")
println("  C₀: [0.1, 5.0]")
println("  b_cap: [0.0, 1.2]")
println("  ε: [0, 365] days")
println()

println("⚠️  CRITICAL ISSUE: ε search range is WRONG!")
println("   Paper's fitted ε = 909 days, but you're only searching up to 365!")
println()

println("Recommended search space:")
println()
println("  OPTION A: If treating C₀ as TOTAL capacity (1.27×10⁵):")
println("    C₀: [5e4, 3e5]")
println("    b_cap: [0.0, 1.2]")
println("    ε: [500, 1500] days  ← CRITICAL: Must include 909!")
println()

println("  OPTION B: If treating C₀ as per-household density:")
println("    C₀: [0.5, 3.0]  (×N_HOUSEHOLDS = 51k to 308k total)")
println("    b_cap: [0.0, 1.2]")
println("    ε: [500, 1500] days  ← CRITICAL: Must include 909!")
println()

# Determine which interpretation is being used
println("To determine which interpretation your code uses:")
println("  Check your entomology.jl get_carrying_capacity() function")
println("  If it multiplies by N_HOUSEHOLDS → use Option B")
println("  If it doesn't multiply → use Option A")
println()

println("\n5. WHAT YOUR FIT RESULTS MEAN:")
println("-"^70)

println("When you get fitted parameters like:")
println("  C₀ = 1.0, b_cap = 0.3, ε = 30.0")
println()
println("This means:")
println("  - Initial carrying capacity: C₀ × N_HOUSEHOLDS = $(round(Int, 1.0 * N_HOUSEHOLDS))")
println("  - Capacity starts growing at day 30 (way too early!)")
println("  - Growth rate: 0.3 per day")
println()
println("But the paper says:")
println("  - Initial capacity: 1.27 × 10⁵")
println("  - Growth starts at day 909 (2.5 years into the data!)")
println("  - Growth rate: 0.3165 per day")
println()

println("⚠️  If your ε fits to small values (<100), it means:")
println("   The optimizer is trying to create early population growth")
println("   because the actual transition (day 909) is outside your search space!")
println()

println("="^70)
println("RECOMMENDATIONS:")
println("="^70)
println()
println("1. UPDATE YOUR run_mosquito_fit.jl IMMEDIATELY:")
println("   Change line 80 to: upper_bounds = [5.0, 1.2, 1500.0]")
println("   Change line 77 to: initial_guess = [1.27, 0.3, 900.0]")
println()
println("2. Check if your entomology.jl multiplies C₀ by N_HOUSEHOLDS")
println("   If YES: use C₀ values around 1-2")
println("   If NO: use C₀ values around 10⁵")
println()
println("3. After refitting, check:")
println("   - Does ε fit to ~900 days?")
println("   - Are predicted MFAI values in range 0.06-0.72?")
println("   - Does the seasonal pattern match?")
println()
println("4. If model still produces too many trapped mosquitoes:")
println("   - Consider reducing α (trap efficiency)")
println("   - Check temperature-dependent mortality rates")
println("   - Verify M(0) = 0.7 × N is appropriate")
println()
println("="^70)

# Create visualization
println("\nGenerating plot...")

p = plot(trap_df.date, mfai_vals,
         marker=:circle,
         markersize=6,
         label="Observed MFAI",
         xlabel="Date",
         ylabel="MFAI (mosquitoes per trap)",
         title="Observed MFAI Over Time",
         legend=:topleft,
         linewidth=2,
         size=(1000, 600))

# Add horizontal line for median
hline!([median(mfai_vals)], label="Median = $(round(median(mfai_vals), digits=3))", 
       linestyle=:dash, color=:red, linewidth=2)

# Add shaded region for paper's fitted ε
if minimum(trap_df.date) isa Date
    start_date = minimum(trap_df.date)
    epsilon_date = start_date + Day(909)
    if epsilon_date <= maximum(trap_df.date)
        vline!([epsilon_date], label="ε = 909 days", 
               linestyle=:dash, color=:green, linewidth=2)
    end
end

savefig(p, "mosquito_diagnostics_plot.png")
println("Plot saved as 'mosquito_diagnostics_plot.png'")

println("\nDiagnostics complete!")
println("\nNext step: Update your run_mosquito_fit.jl with the corrected bounds!")
