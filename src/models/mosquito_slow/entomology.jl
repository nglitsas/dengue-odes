module Entomology

using DataInterpolations
using ..Constants

export get_carrying_capacity, get_rates

# --- 1. The Function to Fit: Carrying Capacity C(t) ---
@inline function smooth_heaviside(x, k)
    z = -2.0 * k * x
    z = clamp(z, -60.0, 60.0)
    return 1.0 / (1.0 + exp(z))
end

function get_carrying_capacity(t, C₀, b_cap, ϵ, t_start)
    # Time since the ramp started (negative before ϵ, positive after)
    dt = t - t_start - ϵ

    # Smooth transition around ϵ (width controlled by k=2.0 → quite sharp)
    u_t = smooth_heaviside(dt, 2.0)

    # Linear growth in capacity after ϵ, flat before
    # b_cap now has units of capacity per time unit (e.g. households per day)
    growth_term = b_cap * dt * u_t   # only grows after transition

    capacity_multiplier = C₀ + growth_term
    final_capacity = capacity_multiplier * Constants.N_HOUSEHOLDS

    return final_capacity
end

# --- 2. Biological Rates (Polynomials) ---
# FIX: Added max(0.0, ...) to all returns to prevent negative rates 
# causing solver explosions at low temperatures.

function oviposition_rate(temp)
    b₀ = -5.3999
    b₁ = 1.800160
    b₂ = -2.12e-1
    b₃ = 1.02e-2
    b₄ = -1.51e-4
    δ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(0.0, δ)  # <--- SAFETY CLAMP
end

function larval_growth(T)
    # Using column sigma_l from your table
    # b0 to b7
     
    b₀ = -1.842702
    b₁ = 8.29e-1
    b₂ = -1.46e-1
    b₃ = 1.31e-2
    b₄ = -6.46e-4
    b₅ = 1.79e-5
    b₆ = -2.62e-7
    b₇ = 1.5e-9
    σl = b₀ + b\1*T + b₂*T^2 + b₃*T^3 + b₄*T^4 + b₅*T^5 + b₆*T^6 + b₇*T^7
    return max(0.01, σl)
end

function pupal_growth(T)
    # Using column sigma_p from your table
    # b0 to b8
    b₀ = 1.9021
    b₁ = -10.31110
    b₂ = 2.050830
    b₃ = -2.24e-1
    b₄ = 1.47e-2
    b₅ = -5.89e-4
    b₆ = 1.41e-5
    b₇ = -1.9e-7
    b₈ = 1.0e-9
    σp = b₀ + b₁*T + b₂*T^2 + b₃*T^3 + b₄*T^4 + b₅*T^5 + b₆*T^6 + b₇*T^7 + b₈*T^8
    return max(0.01, σp)
end

function aquatic_transition(temp)
    σl = larval_growth(temp)
    σp = pupal_growth(temp)
    
    # Combined rate: 1 / (Days_larval + Days_pupal)
    γ = 1.0 / ((1.0/σl) + (1.0/σp))
    return max(0.0, γ)  # <--- SAFETY CLAMP
    # return 0.1
end

function aquatic_mortality(temp)
    b₀ = 2.130
    b₁ = -3.797e-1
    b₂ = 2.457e-2
    b₃ = -6.778e-4
    b₄ = 6.794e-6
    μ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    # Mortality cannot be negative, but we also don't want it to be 0 
    # (nothing lives forever). 
    return max(1e-6, μ) # <--- SAFETY CLAMP
end

function adult_mortality(temp)
    b₀ = 8.69e-1
    b₁ = -1.59e-1
    b₂ = 1.12e-2
    b₃ = -3.41e-4
    b₄ = 3.81e-6
    μₘ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(1e-6, μₘ) # <--- SAFETY CLAMP
end

# --- 3. Rate Accessor ---
function get_rates(t, temp_interp)
    T = temp_interp(t)
    if isnan(T)
        error("Temperature Interpolator returned NaN at t=$t. 
            Check your weather data range!")
    end
    δ = oviposition_rate(T)
    γ_m = aquatic_transition(T)
    μ_a = aquatic_mortality(T)
    μ_m = adult_mortality(T)
    return δ, γ_m, μ_a, μ_m
end

end # module