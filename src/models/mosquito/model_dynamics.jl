#src/models/mosquito/model_dynamics.jl
"""
This file defines the dynamical mosquito capture ODE model that is used to 
estimate highly sensitive mosquito parameters by minimizing the loss between
the model output and observed mosquito trap data 
(data/processed/mosquito_trap-2017_2022.csv). All other parameters are fixed; 
constants are imported from ../constants.jl. and temperature dependent forcing
rates are imported from ../entomology.jl.
"""
module MosquitoModelDynamics

using StaticArrays
using ..Constants
using ..Entomology

export CaptureModel_Fast, CaptureModel, MosquitoModelParams
export compute_mfai_theo, default_u0

# --- 1. Define the Parameter Struct ---
# This ensures type stability and compatibility with your fitter
struct MosquitoModelParams{F}
    C₀::Float64
    bₖ::Float64
    ϵ::Float64
    t_start::Float64
    temp_interp::F
end

# --- 2. The Dynamics Function ---
function CaptureModel_Fast(u, p, t)
    A, M, T = u
    
    # 1. Get Biological Rates
    δₜ, γₘₜ, μₐₜ, μₘₜ = Entomology.get_rates(t, p.temp_interp)

    # 2. Carrying Capacity
    # CRITICAL: We assume Entomology.get_carrying_capacity ALREADY multiplies by N_HOUSEHOLDS.
    # We use the result directly as the Total Capacity.
    C_total = Entomology.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)
    
    # Safety floor
    C_total = max(C_total, 1e-6)

    # 3. Trapping Physics
    # Rate = Alpha * (Traps / Households)
    trap_rate = Constants.ALPHA * (Constants.N_TRAPS / Constants.N_HOUSEHOLDS)
    
    trapping_flow = trap_rate * M
    emergence = γₘₜ * A
    
    # Logistic factor
    logistic_factor = 1.0 - (A / C_total)
    oviposition = Constants.K * δₜ * max(0.0, logistic_factor) * M

    # 4. Derivatives
    dA = oviposition - emergence - (μₐₜ * A)
    dM = emergence - (μₘₜ * M) - trapping_flow
    dT = trapping_flow  # Accumulates total mosquitoes caught

    return SVector(dA, dM, dT)
end

# Alias for compatibility (Non-mutating)
const CaptureModel = CaptureModel_Fast

# --- 3. The Initial Conditions ---
function default_u0(p::MosquitoModelParams)
    # We assume p.C₀ is the DENSITY (Mosquitoes per Household).
    # Therefore, we must scale up to Total Population for the initial state.
    
    # Note: If Entomology already multiplied C₀ by N inside get_carrying_capacity,
    # we still need to do it here manually for the start state because we aren't calling that function yet.
    total_capacity = p.C₀ * Constants.N_HOUSEHOLDS
    
    A_init = 0.85 * total_capacity
    M_init = 0.50 * total_capacity
    T_init = 0
    return SVector{3}(A_init, M_init, T_init)
end

# --- 4. The Corrected MFAI Calculator ---
# This matches the signature your fitter expects and returns a VECTOR.
function compute_mfai_theo(sol, t_steps)
    # Extract the Cumulative Trapped (T) column from the solution at the requested t_steps
    # sol(t, idxs=3) is the interpolation of the Trapped variable
    T_vals = [sol(t)[3] for t in t_steps]
    
    mfai = Vector{Float64}(undef, length(t_steps))
    
    # Handle the first point (assume 0 or steady state, usually 0 catch if t=0)
    mfai[1] = 0.0 
    
    # Loop over the rest
    for i in 2:length(t_steps)
        # 1. How many mosquitoes were caught in this interval?
        delta_catch = T_vals[i] - T_vals[i-1]
        
        # 2. How long was the interval?
        dt = t_steps[i] - t_steps[i-1]
        
        # Prevent division by zero
        dt = max(dt, 1e-5)
        
        # 3. Calculate MFAI: (Catch / Traps) / Days
        # This gives "Mosquitoes per Trap per Day"
        mfai[i] = delta_catch / (Constants.N_TRAPS * dt)
    end
    
    return mfai
end

end # module