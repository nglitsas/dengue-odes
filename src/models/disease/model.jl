module DengueModel

using DifferentialEquations
import DataInterpolations
# FIX: Define the alias explicitly so the struct can find it
const AbstractInterpolation = DataInterpolations.AbstractInterpolation

using ...Shared.TimeUtil
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
    C₀ = 1.33      
    bₖ = 0.3165     
    ϵ  = 989.0          

    # --- BRIERE / THERMAL PARAMETERS ---
    # REQUIRED: EpidemEnto functions (biting_rate, etc.) access these specific fields.
    #biting rate parameters
    b_c::Float64  = 5.22e-4
    b_T0::Float64 = 13.35
    b_Tm::Float64 = 40.08
    
    βₘ_c::Float64  = 4.203e-4
    βₘ_T0::Float64 = 12.22
    βₘ_Tm::Float64 = 37.46

    βₕ_c::Float64  = 1.05e-3
    βₕ_T0::Float64 = 17.05
    βₕ_Tm::Float64 = 35.83

    θₘ_T0::Float64 = 135.0
    θₘ_Tm::Float64 = 14.0

    # --- INITIALIZATION ---
    ϕ  = 0.24              

    # --- DRIVERS ---
    t_start::Float64= 0.0
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
    # C_initial = EpidemEnto.get_carrying_capacity(t0, p.C₀, p.bₖ, p.ϵ)
    N_val = p.N

    C_initial = 0.7*N_val*p.C₀

    M_a0 = 0.35 * C_initial  # 35% aquatic
    M_s0 = 0.50 * C_initial  # 50% susceptible adults
    M_e0 = 43.0
    M_i0 = 43.0

    H_i0 = 62.0
    H_e0 = 62.0
    H_r0 = (1 - p.ϕ)*N_val
    H_s0 = N_val - H_e0 - H_i0 - H_r0
    
    # println("--- Initial State Debug ---")
    # println("Start Time (t0): ", t0)
    # println("Carrying Capacity: ", C_initial)

    return Float64[M_a0, M_s0, M_e0, M_i0, H_s0, H_e0, H_i0, H_r0]
end

# ==========================================
# 5. DYNAMICS
# ==========================================
function DengueModel!(du, u, p::ModelParams, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u

    Mₐ = max(0.0, Mₐ)
    Mₛ = max(0.0, Mₛ)
    Mₑ = max(0.0, Mₑ)
    Mᵢ = max(0.0, Mᵢ)
    Hₛ = max(0.0, Hₛ)
    Hₑ = max(0.0, Hₑ)
    Hᵢ = max(0.0, Hᵢ)
    Hᵣ = max(0.0, Hᵣ)

    # if any(x -> x < 0, u)
    #     date_now = TimeUtil.t_to_date(t)
    #     @warn "NEGATIVE POPULATION at t=$t ($date_now)"
    #     println("  Ma: $Mₐ, Ms: $Mₛ, Me: $Mₑ, Mi: $Mᵢ")
    #     println("  Hs: $Hₛ, He: $Hₑ, Hi: $Hᵢ, Hr: $Hᵣ")
    # end

    # 1. RATES
    # Pass 'p' so Brière functions can find b_c, b_T0, etc.
    r = EpidemEnto.get_rates(t, p.temp_interp, p) 

    # if rand() < 0.001  # Print ~0.1% of timesteps
    #     date_now = TimeUtil.t_to_date(t)
    #     println("t=$t ($date_now): bₜ=$(r.bₜ), βₘₜ=$(r.βₘₜ), βₕₜ=$(r.βₕₜ), θₘₜ=$(r.θₘₜ)")
    # end
    
    # 2. CAPACITY
    C_val = EpidemEnto.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)


    # # ✅ ADD THIS DIAGNOSTIC HERE (after rates and capacity are calculated)
    # if t < p.t_start + 1.0  # Only print at the very start
    #     println("\n=== EQUILIBRIUM CHECK at t=$t ===")
    #     println("Temperature: $(p.temp_interp(t))°C")
    #     println("Capacity: $C_val")
    #     println("Current Mosquitoes: Ma=$Mₐ, Ms=$Mₛ, Total=$(Mₐ + Mₛ + Mₑ + Mᵢ)")
        
    #     M_total = Mₛ + Mₑ + Mᵢ
    #     logistic = (1 - Mₐ / max(C_val, 1e-5))
        
    #     birth_rate = p.k * r.δₜ * logistic * M_total
    #     aquatic_loss = (r.γₘₜ + r.μₐₜ + p.cₐ) * Mₐ
        
    #     println("\nAQUATIC BALANCE:")
    #     println("  Births:           $birth_rate per day")
    #     println("  Losses (death+emergence): $aquatic_loss per day")
    #     println("  Net aquatic:      $(birth_rate - aquatic_loss)")
    #     println("  Ratio:            $(birth_rate / max(aquatic_loss, 1e-10))")
        
    #     adult_recruitment = r.γₘₜ * Mₐ
    #     adult_death = (r.μₘₜ + p.cₘ) * Mₛ
    #     adult_infection = r.bₜ * r.βₕₜ * (Hᵢ / max(p.N, 1.0)) * Mₛ
        
    #     println("\nADULT BALANCE:")
    #     println("  Recruitment:      $adult_recruitment per day")
    #     println("  Deaths:           $adult_death per day")
    #     println("  Getting infected: $adult_infection per day")
    #     println("  Net adults:       $(adult_recruitment - adult_death - adult_infection)")
        
    #     println("\nRATES:")
    #     println("  Oviposition (δ):  $(r.δₜ)")
    #     println("  Emergence (γₘ):   $(r.γₘₜ)")
    #     println("  Aquatic mort (μₐ): $(r.μₐₜ)")
    #     println("  Adult mort (μₘ):  $(r.μₘₜ)")
    #     println("  Biting rate (b):  $(r.bₜ)")
    #     println("="^50)
    # end
    # ✅ END OF DIAGNOSTIC

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