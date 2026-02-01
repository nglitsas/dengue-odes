module DengueODES

# 1. Load Shared Utilities
# Changed from "core/Core.jl" to "shared/Shared.jl"
include("shared/Shared.jl") 

# 2. Load Model Modules
include("models/mosquito/MosquitoCapture.jl")
# include("models/disease/DengueTransmission.jl")

# 3. Export Modules
using .Shared
using .MosquitoCapture
# using .DengueTransmission

export Shared, MosquitoCapture
# , DengueTransmission

end # module
