#src2/mosquito_model.jl
# This file defines the mosquito capture model that is used to estimate mosquito parameters that are then substituted back into the main model

Base.@kwdef struct MosquitoModelParams
    Nⱼ::Float64      # number of mosquito traps (not used in ODE below yet) X
    H₀::Float64      # number of households (not used in ODE below yet) X
    j::Float64       # trap capture rate (not used in ODE below yet) X

    k::Float64       # fraction hatchlings female X
    δₜ::Float64      # oviposition rate X
    μₘₜ::Float64     # adult mosquito mortality X

    μₐₜ::Float64     # aquatic mortality X
    γₘₜ::Float64     # aquatic transition rate (aquatic -> adult susceptible) X

    C::Float64       # mosquito carrying capacity X
    C₀::Float64      # initial carrying capacity (unused in RHS below)
    bₖ::Float64      # carrying capacity coeff (unused in RHS below)
    ϵ::Float64       # carrying capacity threshold (unused in RHS below)
    
end

# State indices (u is length 3)
const IX_A  = 1
const IX_M  = 2
const IX_T  = 3

"""
Default initial condition vector u0 = [A_0, M_0, T_0].

Can be overriden via keyword args.
"""

function default_u0(p::MosquitoModelParams;
    A_0::Real = 0.2p.C,
    M_0::Real = 0.3p.C,
    T_0::Real = 0.0,
   
)
    
    return Float64[
        A_0, M_0, T_0 
    ]
end

function CaptureModel!(du, u, p::MosquitoModelParams, t)
    A, M, T = u

    # transmission rate
    

    # mosquitoes
    du[IX_A] = k*δₜ*(1 - A/C)*M - γₘₜ*A - μₐₜ*A
    du[IX_M] = γₘₜ*A - μₘₜ*M - T
    du[IX_T] = j*(Nⱼ/H₀)*M
    
    return nothing
end