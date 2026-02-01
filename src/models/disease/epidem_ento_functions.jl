module EpidemEnto

using DataInterpolations

# Export the main accessors
export get_carrying_capacity, get_rates

# ==========================================
# 1. HELPER MATH FUNCTIONS
# ==========================================

# Used for Carrying Capacity
@inline function smooth_heaviside(x, k)
    return 1.0 / (1.0 + exp(-2.0 * k * x))
end 

# Generic Briere definition (Type 1) often used for biting/EIP
# Rate = c * T * (T - T0) * sqrt(Tm - T)
@inline function briere_type(t, c, T0, Tm)
    if t < T0 || t > Tm
        return 0.0
    else
        return c * t * (t - T0) * sqrt(max(0.0, Tm - t))
    end
end

# ==========================================
# 2. CARRYING CAPACITY C(t)
# ==========================================
function get_carrying_capacity(t, C₀, bₖ, ϵ)
    dt = t - ϵ
    ramp_active = smooth_heaviside(dt, 2.0)
    return C₀ + bₖ * ramp_active
end

# ==========================================
# 3. BIOLOGICAL RATES (Polynomials)
# ==========================================
# Based on your provided coefficients.
# "Safety Clamps" ensure valid ODE solver steps.

function oviposition_rate(temp)
    b₀ = -5.3999
    b₁ = 1.800160
    b₂ = -2.12e-1
    b₃ = 1.02e-2
    b₄ = -1.51e-4
    δ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(0.0, δ) 
end

function aquatic_transition(temp)
    b₀ = 1.310e-1
    b₁ = -5.723e-1
    b₂ = 1.164e-2
    b₃ = -1.341e-3
    b₄ = 8.723e-5
    γ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(0.0, γ) 
end

function aquatic_mortality(temp)
    b₀ = 2.130
    b₁ = -3.797e-1
    b₂ = 2.457e-2
    b₃ = -6.778e-4
    b₄ = 6.794e-6
    μ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(1e-6, μ) 
end

function adult_mortality(temp)
    b₀ = 8.69e-1
    b₁ = -1.59e-1
    b₂ = 1.12e-2
    b₃ = -3.41e-4
    b₄ = 3.81e-6
    μₘ = b₀ + b₁*temp + b₂*temp^2 + b₃*temp^3 + b₄*temp^4
    return max(1e-6, μₘ) 
end

# ==========================================
# 4. PLACEHOLDERS (To be filled)
# ==========================================

"""
Extrinsic Incubation Rate (θₘ) - Equation 4
Based on a Degree-Day formulation with hard clamps.
"""
function extrinsic_incubation_rate(temp)
    Tm = 14.0   # Threshold temp 
    Ts = 135.0  # Thermal sum 

    # 1. Apply Hard Clamps defined in Methods 
    if temp <= 16.4
        return 0.02
    elseif temp >= 38.0
        return 0.2
    end

    # 2. Calculate Degree-Day Rate
    # Form: (T - Tm) / Ts
    theta = (temp - Tm) / Ts

    # 3. Safety Check (Redundant due to clamps, but good practice)
    return max(0.0, theta)
end

"""
Biting Rate (b)
Often modeled using a Briere function.
"""
function biting_rate(temp)
    # TODO: INSERT BRIERE COEFFICIENTS HERE
    c  = 5.0366e-4
    T0 = 13.35
    Tm = 40.08
    b = briere_type(temp, c, T0, Tm)
    
    return b
    # already clamped in helper function
end

"""
Transmission Probability: Mosquito -> Human (βₘ)
"""
function trans_mosq_human(temp)
    # TODO: INSERT BRIERE COEFFICIENTS HERE
    c  = 6.5093e-4
    T0 = 12.22
    Tm = 37.46
    βₘ = briere_type(temp, c, T0, Tm)
    
    return βₘ
    # already clamped in helper function
end

"""
Transmission Probability: Human -> Mosquito (βₕ)
"""
function trans_human_mosq(temp)
    # TODO: INSERT BRIERE COEFFICIENTS HERE
    c  = 1.0546e-3
    T0 = 17.05
    Tm = 35.83
    βₕ = briere_type(temp, c, T0, Tm)
    
    return βₕ
    # already clamped in helper function
end

# ==========================================
# 5. MAIN ACCESSOR
# ==========================================
"""
get_rates(t, temp_interp)

Queries the temperature at time `t` and returns a NamedTuple
containing all temperature-dependent biological rates.
"""
function get_rates(t, temp_interp)
    T = temp_interp(t)
    
    # 1. Existing Polynomials
    δ   = oviposition_rate(T)
    γₘ  = aquatic_transition(T)
    μₐ  = aquatic_mortality(T)
    μₘ  = adult_mortality(T)

    # 2. New Functional Forms
    θₘ  = extrinsic_incubation_rate(T)
    b   = biting_rate(T)
    βₘ  = trans_mosq_human(T)
    βₕ  = trans_human_mosq(T)

    # Return as a NamedTuple for easy unpacking in the solver
    return (
        δₜ  = δ,
        γₘₜ = γₘ,
        μₐₜ = μₐ,
        μₘₜ = μₘ,
        θₘₜ = θₘ,
        bₜ  = b,
        βₘₜ = βₘ,
        βₕₜ = βₕ
    )
end

end # module