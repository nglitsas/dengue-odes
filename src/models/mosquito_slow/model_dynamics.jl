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

using DifferentialEquations
using ForwardDiff
using ..Constants
using ..Entomology

export MosquitoModelParams, build_params, default_u0, CaptureModel!, IX_A, IX_M, IX_T

# State indices (u is length 3)
const IX_A  = 1 # Aquatic
const IX_M  = 2 # Adult Females
const IX_T  = 3 # Cumulative Trapped

"""
    MosquitoModelParams{I}

Parameters for the mosquito capture model.
Time-varying rates are calculated dynamically using `temp_interp`.
"""
Base.@kwdef struct MosquitoModelParams{I}
    # --- Fixed / Environmental Constants ---
    Nⱼ = Constants.N_TRAPS      # Number of mosquito traps (N_tr)
    H₀ =  Constants.N_HOUSEHOLDS   # Number of households (H_o)
    j = Constants.ALPHA     # Trap capture efficiency (alpha)
    k = Constants.K     # Fraction hatchlings female

    # --- Fitting Parameters for Carrying Capacity C(t) ---
    C₀ = 1.33 # Initial carrying capacity (from paper) per household
    bₖ =  0.3165      # Growth rate (b_cap) (from paper)
    ϵ = 909.0      # Time threshold (epsilon) (from paper)
    t_start::Float64 

    # --- Forcing Data ---
    # Holds the interpolation object (e.g., LinearInterpolation)
    # This allows us to get Temperature T(t) inside the solver
    temp_interp::I   
end

"""
    build_params(temp_interp; fitted_params=nothing)

Create MosquitoModelParams using defaults, optionally overriding [C₀, bₖ, ϵ].
"""
function build_params(temp_interp; t_start::Float64,
    fitted_params::Union{Nothing,AbstractVector{<:Real}}=nothing)
    if fitted_params === nothing
        println("No fitted parameters provided. Using defaults.")
        return MosquitoModelParams(t_start=t_start, temp_interp=temp_interp)  # uses defaults from @kwdef
    end
    @assert length(fitted_params) == 3 "Expected fitted_params = [C₀, bₖ, ϵ]"
    C0, bk, eps = fitted_params
    return MosquitoModelParams(
        C₀ = Float64(C0),
        bₖ = Float64(bk),
        ϵ  = Float64(eps),
        t_start = t_start,
        temp_interp = temp_interp
    )
end

"""
Default initial condition vector u0 = [A, M, T].
"""
function default_u0(p::MosquitoModelParams;
    A₀::Real = 0.85 * p.C₀ * p.H₀, # From paper results
    M₀::Real = 0.7 * Constants.POPULATION, # Scaled to city population
    T₀::Real = 0.0
)
    return Float64[A₀, M₀, T₀]
end

"""
    CaptureModel!(du, u, p, t)

In-place ODE definition.
Dynamically calculates biological rates based on temperature at time t.
"""
function CaptureModel!(du, u, p::MosquitoModelParams, t)
    A, M, T = u
    # if t < p.t_start + 0.00001
    #     @info "Step Check" t Aquatic=A Adults=M
    # end
    #  # --- Input sanity (BEFORE using A,M,T) ---
    # @assert all(x -> isfinite(ForwardDiff.value(x)), u) "u not finite at t=$t"

    # δₜ, γₘₜ, μₐₜ, μₘₜ = Entomology.get_rates(t, p.temp_interp)

    # @assert isfinite(ForwardDiff.value(δₜ))  "δₜ not finite at t=$t"
    # @assert isfinite(ForwardDiff.value(γₘₜ)) "γₘₜ not finite at t=$t"
    # @assert isfinite(ForwardDiff.value(μₐₜ)) "μₐₜ not finite at t=$t"
    # @assert isfinite(ForwardDiff.value(μₘₜ)) "μₘₜ not finite at t=$t"

    # --- 1. Get Dynamic Rates ---
    
    # Calculate all biological rates for current temp T(t)
    # Returns: (oviposition, aquatic_transition, aquatic_mortality, adult_mortality)
    δₜ, γₘₜ, μₐₜ, μₘₜ = Entomology.get_rates(t, p.temp_interp)

    # # 2. EMERGENCY PRINT
    # if t < p.t_start + 0.0001
    #     println("--- ENTOMOLOGY CHECK ---")
    #     println("  Temp: $(p.temp_interp(t))")
    #     println("  Oviposition (δ): $δₜ")
    #     println("  Transition (γ): $γₘₜ")
    #     println("  Larval Death (μ_a): $μₐₜ")
    #     println("  Adult Death (μ_m): $μₘₜ")
    #     println("------------------------")
    # end

    # Calculate Carrying Capacity C(t)
    C_t = Entomology.get_carrying_capacity(t, p.C₀, p.bₖ, p.ϵ, p.t_start)
    
    # Safety clamp for C(t)
    C_t = max(C_t, 1e-6)

    # --- 2. Flows ---

    # Trapping flow: removes from M, adds to T
    # Rate = j * (Traps / Households) * M
    trapping_flow = p.j * (p.Nⱼ / p.H₀) * M

    # Emergence: Aquatic -> Adult
    emergence = γₘₜ * A

    # Oviposition (Births)
    # Logistic term: 1 - A/C(t)
    logistic_factor = 1.0 - (A / C_t)
    
    # Ensure no births if A > C (standard logistic constraint)
    effective_birth_rate = p.k * δₜ * max(0.0, logistic_factor)
    oviposition = effective_birth_rate * M

    # --- 3. Equations ---
    # dA/dt = Births - Emergence - Death
    du[IX_A] = oviposition - emergence - (μₐₜ * A)
    # dM/dt = Emergence - Death - Trapped
    du[IX_M] = emergence - (μₘₜ * M) - trapping_flow
    # dT/dt = Accumulation of trapped mosquitoes
    du[IX_T] = trapping_flow

    # # Check for NaNs
    # if any(x -> !isfinite(ForwardDiff.value(x)), du)
    #     @error "Explosion detected!" t u du C_t δₜ γₘₜ μₐₜ μₘₜ
    #     error("Solver halted due to math explosion")
    # end

    return nothing
end


"""
    compute_mfai(sol, t_points, trap_idx, N_traps)

Converts a cumulative trapped state T(t) from an ODE solution into 
discrete MFAI values (Mosquitoes per Trap) over the intervals defined by t_points.

Arguments:
- `sol`: The DifferentialEquations solution object.
- `t_points`: Vector of times to evaluate (e.g., data collection days).
- `trap_idx`: The index of the Trapped state in the solution vector (usually 3).
- `N_traps`: The number of traps used to normalize the count.
"""
function compute_mfai_theo(sol, t_points, trap_idx, N_traps)
    mfai_values = zeros(length(t_points))
    
    # We assume the cumulative count starts at 0 at t=0.
    # We need to track the previous count to get the delta.
    prev_cum_trapped = 0.0 
    
    # If the first data point isn't t=0, we need to know T at the start of that interval.
    # However, for this specific dataset, we usually just want the delta from 
    # the *previous data point* in the list.
    
    for i in eachindex(t_points)
        t_curr = t_points[i]
        
        # Get cumulative trapped at current time
        current_cum = sol(t_curr)[trap_idx]
        
        # Calculate Delta (New captures since last point)
        newly_trapped = current_cum - prev_cum_trapped
        
        # Normalize
        mfai_values[i] = newly_trapped / N_traps
        
        # Update previous for next iteration
        prev_cum_trapped = current_cum
    end
    
    return mfai_values
end


end # module

