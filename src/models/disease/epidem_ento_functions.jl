module EpidemEnto

using DataInterpolations
# using ..DengueModel.ModelParams

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

function extrinsic_incubation_rate(temp, p)
    if temp <= p.θₘ_T0 || temp >= p.θₘ_Tm
        return 0.125  # median val of thetah from paper
    end
    theta = (temp - p.θₘ_Tm) / (p.θₘ_T0)  # or your formula
    return theta
end

function biting_rate(temp, p)
    return briere_type(temp, p.b_c, p.b_T0, p.b_Tm)
end

function trans_mosq_human(temp, p)
    return briere_type(temp, p.βₘ_c, p.βₘ_T0, p.βₘ_Tm)
end

function trans_human_mosq(temp, p)
    return briere_type(temp, p.βₕ_c, p.βₕ_T0, p.βₕ_Tm)
end

# ==========================================
# 2. MAIN RATE ACCESSOR (combines mosquito + disease rates)
# ==========================================
function get_rates(t, temp_interp, p)
    T = temp_interp(t)
    
    if isnan(T)
        error("Temperature interpolator returned NaN at t=$t")
    end

    # Mosquito rates (imported from Entomology)
    δₜ   = oviposition_rate(T)
    γₘₜ  = aquatic_transition(T)
    μₐₜ  = aquatic_mortality(T)
    μₘₜ  = adult_mortality(T)

    # Disease/transmission rates (now use p)
    θₘₜ  = extrinsic_incubation_rate(T, p)
    bₜ   = biting_rate(T, p)
    βₘₜ  = trans_mosq_human(T, p)
    βₕₜ  = trans_human_mosq(T, p)

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