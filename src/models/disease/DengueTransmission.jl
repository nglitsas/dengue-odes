# DengueTransmission
# This file defines the master module for the full dengue transmission model

module Dengue

using DifferentialEquations
using CSV, DataFrames
using Optim 
using DataInterpolations
using Dates

# 1. Include Shared Modules (Up one level, then into shared)
include(joinpath(@__DIR__, "../shared/Temperature.jl")) 

# 2. Include Local Modules (Same directory)
include("EpidemEnto.jl") 

# 3. Use the Modules
using .Temperature
using .EpidemEnto

# 4. Include Core Logic
include("model.jl")
include("data.jl")
include("Fit.jl")
include("Report.jl")
include("Simulate.jl")

export 
    ModelParams,
    DengueModel!,
    default_u0,
    make_problem,
    simulate,
    load_data,
    fit_lhs,
    plot_fit,
    Temperature,
    EpidemEnto

end