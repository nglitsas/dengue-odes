module DengueModel

using DifferentialEquations
import DataInterpolations
# FIX: Define the alias explicitly so the struct can find it
const AbstractInterpolation = DataInterpolations.AbstractInterpolation

using ..EpidemEnto

export ModelParams, DengueModel!, default_u0, resolve
export IX_Mₐ, IX_Mₛ, IX_Mₑ, IX_Mᵢ, IX_Hₛ, IX_Hₑ, IX_Hᵢ, IX_Hᵣ

# ==========================================
# 1. HELPER
# ==========================================
resolve(x::Function, t) = x(t)
resolve(x, t) = x

const N_HOUSEHOLDS = 4712.0 

# ==========================================
# 2. PARAMETERS
# ==========================================
Base.@kwdef struct ModelParams
    # --- FIXED CONSTANTS (Table 1) ---
    N   = 256088.0         # Human Population
    μₕ  = 3.605e-5         # Human Mortality
    θₕ  = 0.125            # Intrinsic Incubation
    γₕ  = 0.125            # Human Recovery
    k   = 0.5              # Sex Ratio
    cₐ  = 0.0              # Control (Aquatic)
    cₘ  = 0.0              # Control (Adult)

    # --- CARRYING CAPACITY ---
    C₀ = 0.1061326 * N_HOUSEHOLDS             
    bₖ = 0.453566          
    ϵ  = 4910.463          

    # --- BRIERE / THERMAL PARAMETERS ---
    # REQUIRED: EpidemEnto functions (biting_rate, etc.) access these specific fields.
    b_c::Float64  = 2.02e-4
    b_T0::Float64 = 13.35
    b_Tm::Float64 = 40.08
    
    βₘ_c::Float64  = 8.49e-4
    βₘ_T0::Float64 = 17.05
    βₘ_Tm::Float64 = 35.83

    βₕ_c::Float64  = 4.91e-4
    βₕ_T0::Float64 = 12.4
    βₕ_Tm::Float64 = 26.1

    θₘ_T0::Float64 = 10.0
    θₘ_Tm::Float64 = 40.0

    # --- INITIALIZATION ---
    ϕ  = 0.14              

    # --- DRIVERS ---
    temp_interp::AbstractInterpolation 
end

# ==========================================
# 3. STATE INDICES
# ==========================================
const IX_Mₐ = 1
const IX_Mₛ = 2
const IX_Mₑ = 3
const IX_Mᵢ = 4
const IX_Hₛ = 5
const IX_Hₑ = 6
const IX_Hᵢ = 7
const IX_Hᵣ = 8

# ==========================================
# 4. INITIALIZATION
# ==========================================
function default_u0(p::ModelParams, t0::Float64)
    C_initial = EpidemEnto.get_carrying_capacity(t0, p.C₀, p.bₖ, p.ϵ)
    N_val = p.N

    M_a0 = 0.2 * C_initial 
    M_s0 = 0.3 * C_initial 
    M_e0 = 0.0
    M_i0 = 0.0

    H_i0 = 1.0 
    H_e0 = 0.0
    H_r0 = 0.0
    H_s0 = N_val - H_e0 - H_i0 - H_r0
    
    println("--- Initial State Debug ---")
    println("Start Time (t0): ", t0)
    println("Carrying Capacity: ", C_initial)

    return Float64[M_a0, M_s0, M_e0, M_i0, H_s0, H_e0, H_i0, H_r0]
end

# ==========================================
# 5. DYNAMICS
# ==========================================
function DengueModel!(du, u, p::ModelParams, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u

    # 1. RATES
    # Pass 'p' so Brière functions can find b_c, b_T0, etc.
    r = EpidemEnto.get_rates(t, p.temp_interp, p) 
    
    # 2. CAPACITY
    C_val = EpidemEnto.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ)

    # 3. TOTALS
    M = Mₛ + Mₑ + Mᵢ
    H = Hₛ + Hₑ + Hᵢ + Hᵣ 

    # 4. INFECTION
    λₘₕ = r.bₜ * r.βₘₜ * (Mᵢ / max(H, 1.0)) 
    λₕₘ = r.bₜ * r.βₕₜ * (Hᵢ / max(H, 1.0))

    # 5. ODES
    # Aquatic
    du[IX_Mₐ] = p.k * r.δₜ * (1 - Mₐ / max(C_val, 1e-5)) * M - (r.γₘₜ + r.μₐₜ + p.cₐ) * Mₐ
    
    # Mosquitoes
    du[IX_Mₛ] = r.γₘₜ * Mₐ - λₕₘ * Mₛ - (r.μₘₜ + p.cₘ) * Mₛ
    du[IX_Mₑ] = λₕₘ * Mₛ - (r.θₘₜ + r.μₘₜ + p.cₘ) * Mₑ
    du[IX_Mᵢ] = r.θₘₜ * Mₑ - (r.μₘₜ + p.cₘ) * Mᵢ

    # Humans
    du[IX_Hₛ] = p.μₕ * (H - Hₛ) - λₘₕ * Hₛ
    du[IX_Hₑ] = λₘₕ * Hₛ - (p.θₕ + p.μₕ) * Hₑ
    du[IX_Hᵢ] = p.θₕ * Hₑ - (p.γₕ + p.μₕ) * Hᵢ
    du[IX_Hᵣ] = p.γₕ * Hᵢ - p.μₕ * Hᵣ

    return nothing
end

end # module