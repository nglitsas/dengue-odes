# src/model.jl

"""
Dengue transmission model with mosquito aquatic + adult stages, and human SEIR.
Note: parameters marked "temperature dependent" are currently treated as constants.
Later you can replace them with functions of t (via forcing/interpolation).
"""
Base.@kwdef struct ModelParams
    N::Float64       # human population (or locale pop)
    Nⱼ::Float64      # number of mosquito traps (not used in ODE below yet)
    H₀::Float64      # number of households (not used in ODE below yet)
    j::Float64       # trap capture rate (not used in ODE below yet)

    μₕ::Float64      # human mortality rate
    θₕ::Float64      # intrinsic incubation rate (H_e -> H_i)
    γₕ::Float64      # human recovery rate (H_i -> H_r)

    k::Float64       # fraction hatchlings female
    cₐ::Float64      # aquatic control effort
    cₘ::Float64      # adult control effort

    δₜ::Float64      # oviposition rate
    μₘₜ::Float64     # adult mosquito mortality
    μₐₜ::Float64     # aquatic mortality
    γₘₜ::Float64     # aquatic transition rate (aquatic -> adult susceptible)
    θₘₜ::Float64     # extrinsic incubation rate (M_e -> M_i)

    bₜ::Float64      # biting rate (daily bites per mosquito)
    βₘₜ::Float64     # effective mosquito contact rate (mosq gets infected from human)
    βₕₜ::Float64     # effective human contact rate (human gets infected from mosq)

    aᵦ::Float64      # coefficient of biting rate function (unused if bₜ constant)
    αₘ::Float64      # coefficient effective contact mosquito (unused if βₘₜ constant)
    αₕ::Float64      # coefficient effective contact human (unused if βₕₜ constant)

    C::Float64       # mosquito carrying capacity
    C₀::Float64      # initial carrying capacity (unused in RHS below)
    bₖ::Float64      # carrying capacity coeff (unused in RHS below)
    ϵ::Float64       # carrying capacity threshold (unused in RHS below)

    ϕ::Float64       # proportion susceptible initially (for H_s0)
end

# State indices (u is length 8)
const IX_Mₐ  = 1
const IX_Mₛ  = 2
const IX_Mₑ  = 3
const IX_Mᵢ  = 4
const IX_Hₛ  = 5
const IX_Hₑ  = 6
const IX_Hᵢ  = 7
const IX_Hᵣ  = 8

"""
Default initial condition vector u0 = [M_a, M_s, M_e, M_i, H_s, H_e, H_i, H_r].

You can override any of these via keyword args.
"""
function default_u0(p::ModelParams;
    M_a0::Real = 0.2p.C,
    M_s0::Real = 0.3p.C,
    M_e0::Real = 0.0,
    M_i0::Real = 0.0,
    H_e0::Real = 0.0,
    H_i0::Real = 1.0,
    H_r0::Real = 0.0
)
    H_s0 = p.ϕ * p.N - H_e0 - H_i0 - H_r0
    return Float64[
        M_a0, M_s0, M_e0, M_i0,
        H_s0, H_e0, H_i0, H_r0
    ]
end

function DengueModel!(du, u, p::ModelParams, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u

    # totals
    M = Mₛ + Mₑ + Mᵢ
    H = Hₛ + Hₑ + Hᵢ + Hᵣ

    # transmission rate
    λₘₕ = p.bₜ * p.βₕₜ * (Mᵢ / max(H, 1e-12))     # mosq -> human
    λₕₘ = p.bₜ * p.βₘₜ * (Hᵢ / max(H, 1e-12))     # human -> mosq

    # mosquitoes
    du[IX_Mₐ] = p.k * p.δₜ * (1 - Mₐ / p.C) * M - p.γₘₜ*Mₐ - p.μₐₜ*Mₐ - p.cₐ*Mₐ
    du[IX_Mₛ] = p.γₘₜ*Mₐ - (λₕₘ*Mₛ)              - p.μₘₜ*Mₛ - p.cₘ*Mₛ
    du[IX_Mₑ] = (λₕₘ*Mₛ) - p.θₘₜ*Mₑ              - p.μₘₜ*Mₑ - p.cₘ*Mₑ
    du[IX_Mᵢ] = p.θₘₜ*Mₑ                           - p.μₘₜ*Mᵢ - p.cₘ*Mᵢ

    # humans
    du[IX_Hₛ] = p.μₕ*(H - Hₛ) - (λₘₕ*Hₛ)
    du[IX_Hₑ] = (λₘₕ*Hₛ) - (p.θₕ + p.μₕ)*Hₑ
    du[IX_Hᵢ] = p.θₕ*Hₑ   - (p.γₕ + p.μₕ)*Hᵢ
    du[IX_Hᵣ] = p.γₕ*Hᵢ   - p.μₕ*Hᵣ

    return nothing
end