module MosquitoCapture

# -------------------------------------------------------------------
# Load submodules
# -------------------------------------------------------------------

include("Data.jl")
include("constants.jl")
include("entomology.jl")
include("model_dynamics.jl")
include("fitting.jl")
include("Simulate.jl")
include("Report.jl")

# -------------------------------------------------------------------
# Bring submodules into scope
# -------------------------------------------------------------------

using .Data
using .Constants
using .Entomology
using .MosquitoModelDynamics
using .Fitting
using .Simulate
using .Report

# -------------------------------------------------------------------
# Public API
# -------------------------------------------------------------------

export
    # Fitting
    fit_mosquito_model,

    # Simulation
    simulate_mosquito,

    # Reporting
    report_mosquito_fit,

    # Core structs (optional, but useful)
    MosquitoModelParams,
    CaptureModel!,
    default_u0,
    IX_A, IX_M, IX_T

end # module
