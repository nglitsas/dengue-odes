module Entomology

using DataInterpolations

export get_carrying_capacity, get_rates

# --- 1. The Function to Fit: Carrying Capacity C(t) ---
@inline function smooth_heaviside(x, k)
    return 1.0 / (1.0 + exp(-2.0 * k * x))
end 

function get_carrying_capacity(t, C₀, b_cap, ϵ)
    dt = t - ϵ
    ramp_active = smooth_heaviside(dt, 2.0)
    return C₀ + b_cap * ramp_active
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

function aquatic_transition(temp)
    b₀ = 1.310e-1
    b₁ = -5.723e-1
    b₂ = 1.164e-2
    b₃ = -1.341e-3
    b₄ = 8.723e-5
    γ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(0.0, γ)  # <--- SAFETY CLAMP
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
    δ = oviposition_rate(T)
    γ_m = aquatic_transition(T)
    μ_a = aquatic_mortality(T)
    μ_m = adult_mortality(T)
    return δ, γ_m, μ_a, μ_m
end

end # module