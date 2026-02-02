# scripts/sobol_analysis.jl
using GlobalSensitivity, QuasiMonteCarlo, OrdinaryDiffEq, StaticArrays, Plots

# 1. Exact Parameters & Ranges from Fig 5 Caption
param_names = [:δ, :γ_m, :μ_a, :μ_m, :C₀, :ϵ, :b_cap]
lower_bounds = [0.0, 0.0, 0.0234, 0.0301, 10.0, 0.0, 0.001]
upper_bounds = [9.0, 0.2, 0.5, 0.109, 100.0, 1000.0, 1.2]

# 2. Sensitivity Interface
function sobol_interface(p_matrix)
    num_samples = size(p_matrix, 2)
    results = zeros(num_samples)
    
    # Paper uses a long enough window to capture epsilon
    tspan = (0.0, 1000.0) 
    u0 = SVector(100.0, 100.0, 0.0) # A, M, Trapped
    
    # Fixed trapping rate alpha = 0.02 [cite: 137]
    # Trapping calculation: (alpha * Ntr / Ho) * M 
    # Using Table 1 values: Ntr=2412, Ho=102751
    trapping_coeff = 0.02 * (2412 / 102751)

    for i in 1:num_samples
        p = p_matrix[:, i]
        params = (δ=p[1], γ_m=p[2], μ_a=p[3], μ_m=p[4], C₀=p[5], ϵ=p[6], b_cap=p[7])
        
        function gsa_ode(u, p, t)
            A, M, T = u
            
            # Heaviside Linear Growth: C(t) = C0 + bcap * max(0, t - epsilon) 
            Ct = p.C₀ + p.b_cap * max(0.0, t - p.ϵ)
            
            # Equation (2) from paper [cite: 538, 541, 542]
            # k = 0.5 (fraction of female)
            dA = 0.5 * p.δ * (1 - A/Ct) * M - (p.γ_m + p.μ_a) * A
            dM = p.γ_m * A - p.μ_m * M - (trapping_coeff * M)
            dT = trapping_coeff * M
            
            return SVector(dA, dM, dT)
        end
        
        prob = ODEProblem(gsa_ode, u0, tspan, params)
        sol = solve(prob, Tsit5(), save_everystep=false, reltol=1e-4)
        
        # Quantity of Interest: Total Trapped mosquitoes [cite: 140]
        results[i] = sol.u[end][3] 
    end
    return results
end

# 3. Execution
println("--- Running Mirror Fig 5 Sobol Analysis (20 Cores) ---")
bounds = [(lower_bounds[i], upper_bounds[i]) for i in 1:length(lower_bounds)]

# Using Sobol with randomization (scrambling) as suggested by the warning [cite: 139]
res = gsa(sobol_interface, Sobol(), bounds; samples=1000)

# 4. Output and Plotting
st_values = vec(collect(res.ST))
s1_values = vec(collect(res.S1))

println("\nTotal Sensitivity Indices (ST):")
for i in 1:length(param_names)
    println("$(param_names[i]): $(round(st_values[i], digits=4))")
end

p = bar(string.(param_names), st_values, 
        title="Sobol Total Index", 
        ylabel="Sensitivity", legend=false, color=:teal)

savefig("mosquito_capture_sobol_analysis.png")
println("\nPlot saved as mosquito_capture_sobol_analysis.png")