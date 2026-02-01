using DifferentialEquations # Assuming you are using this package

# ==========================================
# 1. THE HELPER
# ==========================================
"""
resolve(param, t)

Returns the value of `param` at time `t`.
- If `param` is a Function, it calls `param(t)`.
- If `param` is a Number, it simply returns it.
"""
resolve(x::Function, t) = x(t)
resolve(x, t) = x

# ==========================================
# 2. THE PARAMETERS (Flexible Types)
# ==========================================
Base.@kwdef struct ModelParams
    # We removed ::Float64 to allow these to be Functions or Numbers
    N            # human population
    Nⱼ = 0.0     # traps (placeholder)
    H₀ = 0.0     # households (placeholder)
    j  = 0.0     # capture rate (placeholder)

    μₕ           # human mortality rate
    θₕ           # intrinsic incubation rate (H_e -> H_i)
    γₕ           # human recovery rate (H_i -> H_r)

    k            # fraction hatchlings female
    cₐ           # aquatic control effort
    cₘ           # adult control effort

    # These are high candidates for time-dependence
    δₜ           # oviposition rate
    μₘₜ          # adult mosquito mortality
    μₐₜ          # aquatic mortality
    γₘₜ          # aquatic transition rate 
    θₘₜ          # extrinsic incubation rate 

    bₜ           # biting rate 
    βₘₜ          # effective contact (mosq -> human)
    βₕₜ          # effective contact (human -> mosq)

    # Coeffs (optional usage)
    aᵦ = 0.0     
    αₘ = 0.0     
    αₕ = 0.0     

    C            # mosquito carrying capacity
    C₀ = 0.0     # initial K (placeholder)
    bₖ = 0.0     # K coeff (placeholder)
    ϵ  = 0.0     # K threshold (placeholder)

    ϕ  = 1.0     # proportion susceptible initially
end

# ==========================================
# 3. STATE INDICES & UTILS
# ==========================================
const IX_Mₐ  = 1
const IX_Mₛ  = 2
const IX_Mₑ  = 3
const IX_Mᵢ  = 4
const IX_Hₛ  = 5
const IX_Hₑ  = 6
const IX_Hᵢ  = 7
const IX_Hᵣ  = 8

function default_u0(p::ModelParams;
    M_a0 = 0.2 * resolve(p.C, 0), # Check C at t=0
    M_s0 = 0.3 * resolve(p.C, 0), 
    M_e0 = 0.0,
    M_i0 = 0.0,
    H_e0 = 0.0,
    H_i0 = 1.0,
    H_r0 = 0.0
)
    # Ensure p.N is a number here for calculation
    N_val = resolve(p.N, 0) 
    H_s0 = p.ϕ * N_val - H_e0 - H_i0 - H_r0
    
    return Float64[
        M_a0, M_s0, M_e0, M_i0,
        H_s0, H_e0, H_i0, H_r0
    ]
end

# ==========================================
# 4. THE ODE MODEL
# ==========================================
function DengueModel!(du, u, p::ModelParams, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u

    # --- A. RESOLVE TIME-DEPENDENT PARAMETERS ---
    # We check if these are functions or numbers at the current time t
    b_val   = resolve(p.bₜ, t)    # Biting rate
    δ_val   = resolve(p.δₜ, t)    # Oviposition
    μ_m_val = resolve(p.μₘₜ, t)   # Adult Mosq Mortality
    θ_m_val = resolve(p.θₘₜ, t)   # Extrinsic Incubation
    C_val   = resolve(p.C, t)     # Carrying Capacity

    # Note: I assumed others are constant for now, but you can 
    # add resolve() calls for any other parameter similarly.

    # --- B. TOTALS ---
    M = Mₛ + Mₑ + Mᵢ
    H = Hₛ + Hₑ + Hᵢ + Hᵣ

    # --- C. TRANSMISSION RATES ---
    # Use b_val (resolved) instead of p.bₜ
    λₘₕ = b_val * p.βₕₜ * (Mᵢ / max(H, 1e-12))     # mosq -> human
    λₕₘ = b_val * p.βₘₜ * (Hᵢ / max(H, 1e-12))     # human -> mosq

    # --- D. DYNAMICS ---
    
    # Aquatic
    # Use δ_val and C_val
    du[IX_Mₐ] = p.k * δ_val * (1 - Mₐ / C_val) * M - p.γₘₜ*Mₐ - p.μₐₜ*Mₐ - p.cₐ*Mₐ
    
    # Adult Mosquitoes
    # Use μ_m_val and θ_m_val
    du[IX_Mₛ] = p.γₘₜ*Mₐ - (λₕₘ*Mₛ)               - μ_m_val*Mₛ - p.cₘ*Mₛ
    du[IX_Mₑ] = (λₕₘ*Mₛ) - θ_m_val*Mₑ             - μ_m_val*Mₑ - p.cₘ*Mₑ
    du[IX_Mᵢ] = θ_m_val*Mₑ                        - μ_m_val*Mᵢ - p.cₘ*Mᵢ

    # Humans
    du[IX_Hₛ] = p.μₕ*(H - Hₛ) - (λₘₕ*Hₛ)
    du[IX_Hₑ] = (λₘₕ*Hₛ) - (p.θₕ + p.μₕ)*Hₑ
    du[IX_Hᵢ] = p.θₕ*Hₑ  - (p.γₕ + p.μₕ)*Hᵢ
    du[IX_Hᵣ] = p.γₕ*Hᵢ  - p.μₕ*Hᵣ

    return nothing
end

# ==========================================
# 5. EXAMPLE USAGE
# ==========================================

# Example: Define a function for Biting Rate (sinusoidal)
function seasonal_biting(t)
    # Oscillates between 0.8 and 1.2
    return 1.0 + 0.2 * sin(2π * t / 365.0)
end

# Define parameters (mixing constants and functions)
params = ModelParams(
    N   = 1000.0,
    μₕ  = 1.0/(70*365),
    θₕ  = 1.0/5.0,
    γₕ  = 1.0/7.0,
    k   = 0.5,
    cₐ  = 0.0,
    cₘ  = 0.0,
    δₜ  = 2.0,           # Constant oviposition
    μₘₜ = 0.1,           # Constant mortality
    μₐₜ = 0.1,
    γₘₜ = 0.1,
    θₘₜ = 0.1,
    bₜ  = seasonal_biting, # <--- PASSING THE FUNCTION HERE
    βₘₜ = 0.3,
    βₕₜ = 0.3,
    C   = 2000.0
)

# Setup Solver
u0 = default_u0(params)
tspan = (0.0, 365.0 * 2) # 2 years
prob = ODEProblem(DengueModel!, u0, tspan, params)

# Solve
sol = solve(prob)

# Plotting (Optional)
# using Plots
# plot(sol, vars=[IX_Hᵢ, IX_Mᵢ], label=["Infected Humans" "Infected Mosq"])