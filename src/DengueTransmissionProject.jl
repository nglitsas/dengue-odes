# CompartmentalModel
# 

module CompartmentalModel

using DifferentialEquations
using CSV, DataFrames
using Optim 

include("model.jl")
include("solver.jl")
include("data.jl")
include("fitting.jl")

export 
    model_params,
    CompartmentalModel!,
    solve_CompartmentalModel,
    load_data,
    fit_CompartmentalModel

end
