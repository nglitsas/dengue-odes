# 1. Define the path to your main module file
# We are in /scripts, so we go up one, then into /src/models/mosquito/
module_path = joinpath(@__DIR__, "..", "src", "models", "mosquito", "MosquitoCapture.jl")

# 2. Load the file and bring the module into scope
include(module_path)
using .MosquitoCapture # Note the dot!

# 3. Rest of your imports
using StaticArrays
using DifferentialEquations
using DataInterpolations

# 1. Mock Data Setup
t_test = collect(1.0:10.0:100.0) # 10 days of data
mfai_obs = rand(length(t_test))  # Random observed data
temp_mock = LinearInterpolation([20.0, 25.0, 22.0], [1.0, 50.0, 100.0])

println("--- Testing Refactored Fitting Pipeline ---")

# 2. Test the new fit_mosquito_model
# We use a very small sample size (10) just to check for errors
#try
result = MosquitoCapture.fit_mosquito_model(
    t_test, 
    mfai_obs, 
    temp_mock; 
    method = :lhs, 
    n_lhs_samples = 20, 
    training_days = 100.0
)
    
    println("\n✅ SUCCESS: The model ran and estimated parameters!")
    println("Estimated C₀: ", result.param[1])
    
#catch e
#    println("\n❌ ERROR DETECTED:")
#    @error e
#    println("\nPossible Culprit: Check if Entomology.get_rates(t, interp) matches.")
#end