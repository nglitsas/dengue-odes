# src/fit.jl
# This file performs the parameter estimation

using Random
using LatinHypercubeSampling
import ..DengueTransmission: ModelParams, simulate, IX_Hₑ

scale(u, lo, hi) = lo + u*(hi-lo)

"Human incidence proxy: new infections entering I_h ≈ θ_h * H_e."
function human_incidence_proxy(sol)
    H_e = Array(sol)[IX_Hₑ, :]
    θ_h = sol.prob.p.θₕ
    return θ_h .* H_e
end

function sse(yhat, y)
    @assert length(yhat) == length(y)
    return sum((yhat .- y).^2)
end

"""
LHS fit over a chosen subset. Returns best_p, best_loss.
`base` is a ModelParams containing fixed values.
`bounds` is Dict for the sampled params.
"""
function fit_lhs(df, base::ModelParams;
                 bounds = Dict(
                   :βₕₜ => (0.01, 1.0),
                   :βₘₜ => (0.01, 1.0),
                   :bₜ  => (0.05, 1.0),
                   :θₘₜ => (0.05, 1.0)
                 ),
                 nsamples::Int=2000,
                 seed::Int=1)

    rng = MersenneTwister(seed)
    t_obs = Float64.(df.t)
    y_obs = Float64.(df.cases)
    tspan = (minimum(t_obs), maximum(t_obs))

    keys_ = collect(keys(bounds))
    d = length(keys_)
    U = sample(rng, LatinHypercube(nsamples, d))

    best_loss = Inf
    best_p = base

    for k in 1:nsamples
        # start from base, update sampled params
        new = Dict{Symbol,Float64}()
        for (i, sym) in enumerate(keys_)
            lo, hi = bounds[sym]
            new[sym] = scale(U[k,i], lo, hi)
        end

        p = ModelParams(;  # keyword constructor thanks to Base.@kwdef
            N=base.N, Nⱼ=base.Nⱼ, H₀=base.H₀, j=base.j,
            μₕ=base.μₕ, θₕ=base.θₕ, γₕ=base.γₕ,
            k=base.k, cₐ=base.cₐ, cₘ=base.cₘ,
            δₜ=base.δₜ, μₘₜ=base.μₘₜ, μₐₜ=base.μₐₜ, γₘₜ=base.γₘₜ,
            θₘₜ=get(new, :θₘₜ, base.θₘₜ),
            bₜ=get(new, :bₜ, base.bₜ),
            βₘₜ=get(new, :βₘₜ, base.βₘₜ),
            βₕₜ=get(new, :βₕₜ, base.βₕₜ),
            aᵦ=base.aᵦ, αₘ=base.αₘ, αₕ=base.αₕ,
            C=base.C, C₀=base.C₀, bₖ=base.bₖ, ϵ=base.ϵ, ϕ=base.ϕ
        )

        sol = simulate(p, tspan; saveat=t_obs)
        yhat = human_incidence_proxy(sol)
        L = sse(yhat, y_obs)

        if L < best_loss
            best_loss = L
            best_p = p
        end
    end

    return best_p, best_loss
end
