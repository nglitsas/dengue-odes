module Entomology

using DataInterpolations

export get_carrying_capacity, get_rates

# --- 1. The Function to Fit: Carrying Capacity C(t) ---
# Equation (6): Smooth Heaviside function
# Sophie note: consider switching this to our own function or an ML function
# Robust, differentiable approximation of Heaviside
# k determines the steepness (higher = sharper step)
@inline function smooth_heaviside(x, k)
    return 1.0 / (1.0 + exp(-2.0 * k * x))
end 
function get_carrying_capacity(t, C₀, b_cap, ϵ)
    dt = t - ϵ
    ramp_active = smooth_heaviside(dt, 2.0)

    # Clamp the (t - epsilon) part with a smooth ReLU so it doesn't contribute 
    # negative values even if the heaviside leaks slightly.
    return C₀ + b_cap * dt * ramp_active
end

# --- 2. Biological Rates (Polynomials) ---
# These are inputs. You need to fill in the polynomial coefficients 
# from Yang et al. [32, 28].

function oviposition_rate(temp)
    b₀ = -5.3999
    b₁ = 1.800160
    b₂ = -2.12e-1
    b₃ = 1.02e-2
    b₄ = -1.51e-4
    δ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return δ  
end

function aquatic_transition(temp)
    # placeholder: gamma_m(T) NICK TO ADD
    b₀ = -5.3999
    b₁ = 1.800160
    b₂ = -2.12e-1
    b₃ = 1.02e-2
    b₄ = -1.51e-4
    δ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return δ  
end

function aquatic_mortality(temp)
    # placeholder: mu_a(T) NICK TO ADD
    b₀ = -5.3999
    b₁ = 1.800160
    b₂ = -2.12e-1
    b₃ = 1.02e-2
    b₄ = -1.51e-4
    δ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return δ  
end

function adult_mortality(temp)
    b₀ = 8.69e-1
    b₁ = -1.59e-1
    b₂ = 1.12e-2
    b₃ = -3.41e-4
    b₄ = 3.81e-6
    μₘ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return μₘ
end

# --- 3. Rate Accessor ---
# Helper to get all rates at once, given a temperature interpolation object
# defined in src/core/Temperature.jl
function get_rates(t, temp_interp)
    # Get temp at this time t using the interpolator 
    T = temp_interp(t)
    # Calculate rates
    δ = oviposition_rate(Tₓ)
    γ_m = aquatic_transition(Tₓ)
    μ_a = aquatic_mortality(Tₓ)
    μ_m = adult_mortality(Tₓ)
    
    return δ, γ_m, μ_a, μ_m
end

end