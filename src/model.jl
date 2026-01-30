# base_model
# This file defines the compartmental model for dengue transmission that accounts for climate effects

struct model_params
    αₘ::Float64
    γₕ::Float64
    σₘ::Float64
    σₕ::Float64
    cₘ::Float64
    cₐ::Float64
    μₘ::Float64
    μₕ::Float64
    λₘₕ::Float64
    λₕₘ::Float64
end

function CompartmentalModel!(du, u, p::model_params, t)
    Mₐ, Mₛ, Mₑ, Mᵢ, Hₛ, Hₑ, Hᵢ, Hᵣ = u 
    αₘ, γₕ, σₘ, σₕ, cₘ, cₐ, μₘ, μₕ, λₘₕ, λₕₘ = αₘ.p, γₕ.p, σₘ.p, σₕ.p, cₘ.p, cₐ.p, μₘ.p, μₕ.p, λₘₕ.p, λₕₘ.p


    du[1] = 
    du[2] = 
    du[3] = 
    du[4] =
    du[5] = 
    du[6] =
    du[7] = 
    du[8] = 