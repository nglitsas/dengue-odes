using DifferentialEquations
import DataInterpolations
const AbstractInterpolation = DataInterpolations.AbstractInterpolation

# ==========================================
# 1. THE HELPER
# ==========================================
resolve(x::Function, t) = x(t)
resolve(x, t) = x

# struct_defs.jl

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
    ϵ  = 909.0             

    # --- INITIALIZATION ---
    ϕ  = 0.14              

    # --- DRIVERS ---
    # We replace the individual rate placeholders with just the temperature data
    temp_interp::AbstractInterpolation 
end

# ==========================================
# 3. STATE INDICES & INITIALIZATION
# ==========================================
const IX_Mₐ = 1
const IX_Mₛ = 2
const IX_Mₑ = 3
const IX_Mᵢ = 4
const IX_Hₛ = 5
const IX_Hₑ = 6
const IX_Hᵢ = 7
const IX_Hᵣ = 8

function default_u0(p::ModelParams, t0::Float64)
    # Use the actual start time t0 instead of 0.0
    # Calculate the initial carrying capacity using the parameters in p
    # We use p.C₀, p.bₖ, and p.ϵ because those ARE in your ModelParams
    C_initial = get_carrying_capacity(t0, p.C₀, p.bₖ, p.ϵ)

    N_val = p.N

    # Note: Paper assumes specific start proportions, usually mostly susceptible
    M_a0 = 0.2 * C_initial  # Arbitrary initialization for aquatic
    M_s0 = 0.3 * C_initial  # Arbitrary initialization for adults
    M_e0 = 0.0
    M_i0 = 0

    # Human Initial states
    # If starting with an outbreak, H_i0 > 0
    H_i0 = 1.0 
    H_e0 = 0.0
    H_r0 = 0.0
    # Remaining population is Susceptible
    H_s0 = N_val - H_e0 - H_i0 - H_r0
    
    # Diagnostic for intial values
    println("--- Initial State Debug ---")
    println("Start Time (t0): ", t0)
    println("Carrying Capacity: ", C_initial)
    println("Initial Humans (Ih): ", 1.0)
    println("Initial Mosquitoes (Ms): ", 0.3 * C_initial)

    return Float64[
        M_a0, M_s0, M_e0, M_i0,
        H_s0, H_e0, H_i0, H_r0
    ]

end

# dynamics.jl
function DengueModel!(du, u, p::ModelParams, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u

    # 1. GET ALL BIOLOGICAL RATES FOR CURRENT TEMP
    # This calls your module
    r = EpidemEnto.get_rates(t, p.temp_interp) 

    # 2. CALCULATE CARRYING CAPACITY
    C_val = EpidemEnto.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ)

    # 3. TOTALS
    M = Mₛ + Mₑ + Mᵢ
    H = Hₛ + Hₑ + Hᵢ + Hᵣ 

    # 4. FORCE OF INFECTION
    # r.bₜ is the Biting Rate from your named tuple
    # r.βₘₜ is Mosq->Human transmission
    # r.βₕₜ is Human->Mosq transmission
    λₘₕ = r.bₜ * r.βₘₜ * (Mᵢ / max(H, 1.0)) 
    λₕₘ = r.bₜ * r.βₕₜ * (Hᵢ / max(H, 1.0))

    # 5. DYNAMICS
    # Aquatic
    # r.δₜ = Oviposition, r.γₘₜ = Transition, r.μₐₜ = Aquatic Mortality
    du[IX_Mₐ] = p.k * r.δₜ * (1 - Mₐ / C_val) * M - (r.γₘₜ + r.μₐₜ + p.cₐ) * Mₐ
    
    # Mosquitoes
    # r.μₘₜ = Adult Mortality, r.θₘₜ = Extrinsic Incubation
    du[IX_Mₛ] = r.γₘₜ * Mₐ - λₕₘ * Mₛ - (r.μₘₜ + p.cₘ) * Mₛ
    du[IX_Mₑ] = λₕₘ * Mₛ - (r.θₘₜ + r.μₘₜ + p.cₘ) * Mₑ
    du[IX_Mᵢ] = r.θₘₜ * Mₑ - (r.μₘₜ + p.cₘ) * Mᵢ

    # Humans (Standard SEIR)
    du[IX_Hₛ] = p.μₕ * (H - Hₛ) - λₘₕ * Hₛ
    du[IX_Hₑ] = λₘₕ * Hₛ - (p.θₕ + p.μₕ) * Hₑ
    du[IX_Hᵢ] = p.θₕ * Hₑ - (p.γₕ + p.μₕ) * Hᵢ
    du[IX_Hᵣ] = p.γₕ * Hᵢ - p.μₕ * Hᵣ

    return nothing
end