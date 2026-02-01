module DengueUDE

using DifferentialEquations
using SciMLSensitivity
using Optimization
using OptimizationOptimisers   # ADAM
using OptimizationOptimJL      # LBFGS
using ComponentArrays
using Lux
using Random
using Printf
using Statistics

using ..Model  # DengueModel!, ModelParams, default_u0, IX_*
using ..EpidemEnto  # get_rates, get_carrying_capacity