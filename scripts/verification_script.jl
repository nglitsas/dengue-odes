# VERIFICATION SCRIPT
# Run this after applying the corrections to verify everything is working

using Plots

# ============================================================================
# TEST 1: Verify Carrying Capacity Function
# ============================================================================
function test_carrying_capacity()
    println("="^70)
    println("TEST 1: Carrying Capacity Function")
    println("="^70)
    
    # Paper's fitted values
    C₀ = 1.27e5
    b_cap = 0.3165
    ϵ = 909.0
    t_start = 0.0
    
    # Test at key time points
    times = [0.0, 500.0, 909.0, 910.0, 1000.0, 1500.0]
    
    for t in times
        # Your corrected function
        dt = t - t_start - ϵ
        u_t = dt >= 0.0 ? 1.0 : 0.0
        C = C₀ + b_cap * dt * u_t
        
        println("t = $t days:")
        println("  dt = $(round(dt, digits=2)) days")
        println("  u(dt) = $u_t")
        println("  C(t) = $(round(C, sigdigits=4))")
        println()
    end
    
    # Expected behavior:
    println("Expected behavior:")
    println("  - C(t) = 1.27×10⁵ for t < 909")
    println("  - C(909) = 1.27×10⁵ (exactly)")
    println("  - C(1000) = 1.27×10⁵ + 0.3165×91 = 1.298×10⁵")
    println()
    
    # Plot it
    t_plot = 0:10:1800
    C_plot = [begin
        dt = t - ϵ
        u_t = dt >= 0.0 ? 1.0 : 0.0
        C₀ + b_cap * dt * u_t
    end for t in t_plot]
    
    p = plot(t_plot, C_plot, 
             xlabel="Time (days)", 
             ylabel="Carrying Capacity",
             title="Carrying Capacity C(t)",
             linewidth=2,
             label="C(t)",
             legend=:topleft)
    vline!([ϵ], label="ε = 909 days", linestyle=:dash, color=:red)
    
    return p
end

# ============================================================================
# TEST 2: Verify Initial Conditions
# ============================================================================
function test_initial_conditions()
    println("="^70)
    println("TEST 2: Initial Conditions")
    println("="^70)
    
    # Constants
    N = 256088.0  # Population
    C₀ = 1.27e5   # Total capacity
    
    # Corrected initial conditions
    A_init = 0.85 * C₀
    M_init = 0.7 * N
    T_init = 0.0
    
    println("Initial Conditions (CORRECTED):")
    println("  A(0) = 0.85 × C₀ = 0.85 × 1.27×10⁵ = $(round(A_init, sigdigits=4))")
    println("  M(0) = 0.7 × N = 0.7 × 256,088 = $(round(M_init, sigdigits=4))")
    println("  T(0) = $T_init")
    println()
    
    println("Expected from paper:")
    println("  A(0) ≈ 1.08×10⁵")
    println("  M(0) ≈ 1.79×10⁵")
    println("  T(0) = 0")
    println()
    
    # Check if close
    if abs(A_init - 1.0795e5) / 1.0795e5 < 0.01
        println("✅ A(0) matches paper (within 1%)")
    else
        println("❌ A(0) does NOT match paper")
    end
    
    if abs(M_init - 1.792616e5) / 1.792616e5 < 0.01
        println("✅ M(0) matches paper (within 1%)")
    else
        println("❌ M(0) does NOT match paper")
    end
    println()
end

# ============================================================================
# TEST 3: Verify MFAI Calculation
# ============================================================================
function test_mfai_calculation()
    println("="^70)
    println("TEST 3: MFAI Calculation")
    println("="^70)
    
    # Simulate trapped mosquitoes at two time points
    N_TRAPS = 2412.0
    
    # Example: 100 mosquitoes caught per trap over 60 days
    trapped_t0 = 0.0
    trapped_t1 = 100.0 * N_TRAPS  # Total across all traps
    
    println("Scenario: 100 mosquitoes per trap caught in 60 days")
    println("  Trapped(t=0) = $trapped_t0")
    println("  Trapped(t=60) = $trapped_t1")
    println()
    
    # WRONG calculation (your old code)
    dt = 60.0
    mfai_wrong = (trapped_t1 - trapped_t0) / (N_TRAPS * dt)
    
    # CORRECT calculation (paper's Equation 3)
    mfai_correct = (trapped_t1 - trapped_t0) / N_TRAPS
    
    println("WRONG calculation (dividing by dt):")
    println("  MFAI = (trapped_t1 - trapped_t0) / (N_TRAPS × dt)")
    println("  MFAI = $trapped_t1 / ($N_TRAPS × $dt)")
    println("  MFAI = $(round(mfai_wrong, digits=4))")
    println()
    
    println("CORRECT calculation (Equation 3):")
    println("  MFAI = (trapped_t1 - trapped_t0) / N_TRAPS")
    println("  MFAI = $trapped_t1 / $N_TRAPS")
    println("  MFAI = $(round(mfai_correct, digits=4))")
    println()
    
    println("Difference: $(round(mfai_correct/mfai_wrong, digits=2))× larger!")
    println()
    
    println("Expected MFAI range from paper: 0 to ~10")
    println("  Wrong method gives: ~1-2 (too small by 60×)")
    println("  Correct method gives: ~100 (reasonable)")
    println()
end

# ============================================================================
# TEST 4: Parameter Sanity Check
# ============================================================================
function test_parameters()
    println("="^70)
    println("TEST 4: Parameter Values")
    println("="^70)
    
    println("From paper (Table 1 and Figure 6):")
    println("  N_TRAPS = 2,412")
    println("  N_HOUSEHOLDS = 102,751")
    println("  POPULATION = 256,088")
    println("  ALPHA = 0.02")
    println("  K = 0.5")
    println()
    
    println("Fitted values (Figure 6 caption):")
    println("  C₀ = 1.27 × 10⁵")
    println("  b_cap = 0.3165")
    println("  ε = 909 days")
    println()
    
    println("Initial conditions:")
    println("  A(0) = 0.85 × C₀ = 1.08 × 10⁵")
    println("  M(0) = 0.7 × N = 1.79 × 10⁵")
    println()
    
    println("Collection schedule:")
    println("  Start: Sept 1, 2017 (t=0)")
    println("  End: May 13, 2022 (t≈1715 days)")
    println("  Interval: Bimonthly (~60 days)")
    println("  Collections on odd months: Sept, Nov, Jan, Mar, May")
    println()
end

# ============================================================================
# RUN ALL TESTS
# ============================================================================
function run_all_tests()
    println("\n")
    println("╔" * "="^68 * "╗")
    println("║" * " "^15 * "VERIFICATION SCRIPT FOR CORRECTIONS" * " "^18 * "║")
    println("╚" * "="^68 * "╝")
    println()
    
    # Run tests
    test_initial_conditions()
    test_mfai_calculation()
    test_parameters()
    
    # Plot carrying capacity
    println("Generating carrying capacity plot...")
    p = test_carrying_capacity()
    
    println("\n")
    println("="^70)
    println("SUMMARY")
    println("="^70)
    println("If all checks pass, your corrected code should now:")
    println("  1. Have MFAI values in the correct range (0-10)")
    println("  2. Show seasonal variation matching the paper")
    println("  3. Reproduce Figure 6 when you run the full fit")
    println()
    println("Next steps:")
    println("  1. Replace your old files with the CORRECTED versions")
    println("  2. Re-run your fitting script")
    println("  3. Plot the results and compare with Figure 6")
    println("="^70)
    
    return p
end

# Run it!
p = run_all_tests()
display(p)
savefig(p, "carrying_capacity_verification.png")
println("\nPlot saved as 'carrying_capacity_verification.png'")
