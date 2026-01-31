# src/simulation.jl
# This file sets up the ode solver and wrappers

using DifferentialEquations
import ..DengueTransmission: DengueModel!, ModelParams, default_u0

function make_problem(p::ModelParams, tspan::Tuple{<:Real,<:Real}; u0=default_u0(p))
    ODEProblem(DengueModel!, u0, tspan, p)
end

function simulate(p::ModelParams, tspan; saveat=1.0, u0=default_u0(p))
    prob = make_problem(p, tspan; u0=u0)
    solve(prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
end