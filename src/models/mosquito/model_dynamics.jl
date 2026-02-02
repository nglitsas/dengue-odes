"""
CORRECTED VERSION of model_dynamics.jl

Key fixes:
1. MFAI calculation: Removed division by dt (Equation 3)
2. Initial conditions: M_init now uses POPULATION, not capacity
"""
module MosquitoModelDynamics

using StaticArrays
using ..Constants
using ..Entomology

export CaptureModel_Fast, CaptureModel, MosquitoModelParams
export compute_mfai_theo, default_u0

# --- 1. Define the Parameter Struct ---
struct MosquitoModelParams{F}
    C₀::Float64
    bₖ::Float64
    ϵ::Float64
    t_start::Float64
    temp_interp::F
end

# --- 2. The Dynamics Function ---
function CaptureModel_Fast(u, p, t)
    A, M, T = u
    
    # 1. Get Biological Rates
    δₜ, γₘₜ, μₐₜ, μₘₜ = Entomology.get_rates(t, p.temp_interp)

    # Numerical safety
    δₜ  = max(0.0, δₜ)
    γₘₜ = max(0.0, γₘₜ)
    μₐₜ = max(0.0, μₐₜ)
    μₘₜ = max(0.0, μₘₜ)

    # 2. Carrying Capacity
    # ✅ FIX: Entomology.get_carrying_capacity should return TOTAL capacity
    # not per-household density. See corrected entomology.jl
    C_total = Entomology.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)
    C_total = max(C_total, 1e-6)

    # 3. Trapping Physics
    trap_rate = Constants.ALPHA * (Constants.N_TRAPS / Constants.N_HOUSEHOLDS)
    trapping_flow = trap_rate * M
    emergence = γₘₜ * A
    
    # Logistic factor
    logistic_factor = 1.0 - (A / C_total)
    oviposition = Constants.K * δₜ * max(0.0, logistic_factor) * M

    # 4. Derivatives
    dA = oviposition - emergence - (μₐₜ * A)
    dM = emergence - (μₘₜ * M) - trapping_flow
    dT = trapping_flow

    return SVector(dA, dM, dT)
end

# Alias for compatibility
const CaptureModel = CaptureModel_Fast

# --- 3. The Initial Conditions ---
function default_u0(p::MosquitoModelParams)
    # C₀ is now TOTAL capacity (not per-household)
    A_init = 0.85 * p.C₀
    M_init = 0.7 * Constants.POPULATION
    T_init = 0.0
    return SVector{3}(A_init, M_init, T_init)
end

# --- 4. The MFAI Calculator ---
function compute_mfai_theo(sol, t_steps)
    """
    ✅ CORRECTED: Equation (3) from paper
    MFAI_theo(tb) = (Trapped(tb) - Trapped(tb-1)) / Ntr
    
    NO DIVISION BY dt! The MFAI is the total catch between collections 
    divided by the number of traps.
    """
    T_vals = [sol(t)[3] for t in t_steps]
    
    mfai = Vector{Float64}(undef, length(t_steps))
    mfai[1] = 0.0 
    
    for i in 2:length(t_steps)
        # Total mosquitoes caught in this interval (across all traps)
        delta_catch = T_vals[i] - T_vals[i-1]
        
        # ✅ CORRECTED: Just divide by number of traps
        # The paper explicitly states this is NOT per day
        mfai[i] = delta_catch / Constants.N_TRAPS
    end
    
    return mfai
end

end # module