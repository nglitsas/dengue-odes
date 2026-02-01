# DengueTransmission
# This file defines the master module for the full dengue transmission model

module Dengue

using DifferentialEquations
using CSV, DataFrames
using Optim 
using DataInterpolations
using Dates
using RollingFunctions

# 1. Access Shared Modules from the Parent (DengueODES -> Shared)
using ..Shared.Temperature
using ..Shared.TimeUtil
using ..MosquitoCapture.Entomology

# 2. Include Local Modules (Same directory)
include("epidem_ento_functions.jl") 

# 3. Include Core Logic
include("model.jl")
include("data.jl")
include("Simulate.jl")
include("Fit.jl")
include("Report.jl")

# 3. Use the Modules
using .Temperature
using .EpidemEnto
using .DengueModel
using .Simulate


export 
    ModelParams,
    DengueModel!,
    default_u0,
    make_problem,
    Simulate,
    load_data,
    fit_lhs,
    plot_fit,
    Temperature,
    EpidemEnto

end