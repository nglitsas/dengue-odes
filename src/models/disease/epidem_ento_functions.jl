module EpidemEnto

using DataInterpolations

# Import ONLY the mosquito-specific rates and carrying capacity from Entomology
using DengueODES.MosquitoCapture.Entomology: get_carrying_capacity,
                    oviposition_rate,
                    aquatic_transition,
                    aquatic_mortality,
                    adult_mortality

# Export everything so dengue code can use EpidemEnto.get_rates(...) cleanly
export get_carrying_capacity, get_rates

# ==========================================
# 1. DISEASE/TRANSMISSION-SPECIFIC FUNCTIONS (Briere)
# ==========================================

@inline function briere_type(t, c, T0, Tm)
    if t < T0 || t > Tm
        return 0.0
    else
        return c * t * (t - T0) * sqrt(max(0.0, Tm - t))
    end
end

function extrinsic_incubation_rate(temp)
    Tm = 14.0    # Threshold temp 
    Ts = 135.0   # Thermal sum 

    if temp <= 16.4
        return 0.02
    elseif temp >= 38.0
        return 0.2
    end

    theta = (temp - Tm) / Ts
    return max(0.0, theta)
end

function biting_rate(temp)
    c  = 5.0366e-4
    T0 = 13.35
    Tm = 40.08
    return briere_type(temp, c, T0, Tm)
end

function trans_mosq_human(temp)
    c  = 6.5093e-4
    T0 = 12.22
    Tm = 37.46
    return briere_type(temp, c, T0, Tm)
end

function trans_human_mosq(temp)
    c  = 1.0546e-3
    T0 = 17.05
    Tm = 35.83
    return briere_type(temp, c, T0, Tm)
end

# ==========================================
# 2. MAIN RATE ACCESSOR (combines mosquito + disease rates)
# ==========================================
function get_rates(t, temp_interp)
    T = temp_interp(t)
    
    if isnan(T)
        error("Temperature interpolator returned NaN at t=$t")
    end

    # Mosquito rates (imported from Entomology)
    δₜ   = oviposition_rate(T)
    γₘₜ  = aquatic_transition(T)
    μₐₜ  = aquatic_mortality(T)
    μₘₜ  = adult_mortality(T)

    # Disease/transmission rates (defined here)
    θₘₜ  = extrinsic_incubation_rate(T)
    bₜ   = biting_rate(T)
    βₘₜ  = trans_mosq_human(T)
    βₕₜ  = trans_human_mosq(T)

    return (
        δₜ   = δₜ,
        γₘₜ  = γₘₜ,
        μₐₜ  = μₐₜ,
        μₘₜ  = μₘₜ,
        θₘₜ  = θₘₜ,
        bₜ   = bₜ,
        βₘₜ  = βₘₜ,
        βₕₜ  = βₕₜ
    )
end

end # module EpidemEnto