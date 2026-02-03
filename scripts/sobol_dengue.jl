# scripts/sobol_dengue_analysis.jl
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GlobalSensitivity, QuasiMonteCarlo, OrdinaryDiffEq, StaticArrays, Plots, Logging
using DengueODES
using DengueODES.Shared.Temperature
using DengueODES.Dengue.DengueModel
using DengueODES.Dengue.EpidemEnto
using Printf
using Base.Threads   # ✅ ADDED

println("="^60)
println("DENGUE TRANSMISSION MODEL: SOBOL SENSITIVITY ANALYSIS")
println("="^60)

# 1. Load Temperature Data (required for all simulations)
println("\nLoading temperature data...")
weather_df = Temperature.get_weather_data()
temp_interp = Temperature.get_temperature_interpolator(weather_df)

# 2. Define Parameters & Ranges
param_names = [:C₀, :b_c, :βₘ_c, :βₕ_c, :ϕ, :θₕ, :γₕ]

# Parameter bounds based on literature
lower_bounds = [
    0.5,      # C₀: carrying capacity (per household)
    0.0005,     # b_c: biting rate Briere coefficient
    0.0005,     # βₘ_c: mosquito transmission coefficient
    0.0005,     # βₕ_c: human transmission coefficient
    0.05,     # ϕ: initial immune fraction (5-40%)
    0.08,     # θₕ: human incubation rate (1/20 to 1/5 days)
    0.05      # γₕ: human recovery rate (1/20 to 1/5 days)
]

upper_bounds = [
    18,      # C₀
    0.0015,     # b_c
    0.0015,     # βₘ_c
    0.0015,     # βₕ_c
    0.40,     # ϕ
    0.33,     # θₕ
    0.25      # γₕ
]

println("\nParameter ranges:")
for i in 1:length(param_names)
    @printf("  %-6s: [%.4f, %.4f]\n", param_names[i], lower_bounds[i], upper_bounds[i])
end

const _DEVNULL = open(devnull, "w")

function quietly(f::Function)
    # Silences:
    # - println/print (stdout)
    # - many warnings and solver output (stderr)
    # - @warn/@info logs via Logging (NullLogger)
    redirect_stdout(_DEVNULL) do
        redirect_stderr(_DEVNULL) do
            with_logger(NullLogger()) do
                return f()
            end
        end
    end
end

# 3. Sensitivity Interface Function
function sobol_interface(p_matrix)
    num_samples = size(p_matrix, 2)
    results = zeros(num_samples)

    # Simulation settings
    t_start = 6250.0  # Start time (2017)
    tspan = (t_start, t_start + 730.0)  # 2 years

    println("\nRunning $num_samples parameter combinations...")
    println("Threads available: $(Threads.nthreads())")  # ✅ ADDED (tiny)

    # ✅ ADDED: thread-safe progress counter
    progress = Threads.Atomic{Int}(0)

    # ✅ CHANGED: threaded loop
    Threads.@threads for i in 1:num_samples
        # progress printing (avoid too much interleaved spam)
        pi = Threads.atomic_add!(progress, 1) + 1
        if pi % 50 == 0
            # printing from multiple threads may interleave, but ok
            @printf("  Progress: %d/%d (%.1f%%)\n", pi, num_samples, 100*pi/num_samples)
        end

        try
            # Extract sampled parameters
            C₀_sample  = p_matrix[1, i]
            b_c_sample = p_matrix[2, i]
            βₘ_c_sample = p_matrix[3, i]
            βₕ_c_sample = p_matrix[4, i]
            ϕ_sample   = p_matrix[5, i]
            θₕ_sample  = p_matrix[6, i]
            γₕ_sample  = p_matrix[7, i]

            # Create parameter struct with sampled values
            params = DengueModel.ModelParams(
                C₀ = C₀_sample,
                b_c = b_c_sample,
                βₘ_c = βₘ_c_sample,
                βₕ_c = βₕ_c_sample,
                ϕ = ϕ_sample,
                θₕ = θₕ_sample,
                γₕ = γₕ_sample,
                t_start = t_start,
                temp_interp = temp_interp
            )

            # Get initial conditions
            u0 = DengueModel.default_u0(params, t_start)

            # Solve ODE
            prob = ODEProblem(DengueModel.DengueModel!, u0, tspan, params)

            sol = quietly() do
                solve(prob, Rosenbrock23(),
                    save_everystep=false,
                    reltol=1e-4,
                    abstol=1e-6,
                    maxiters=1e6)
            end

            if sol.retcode == ReturnCode.Success
                # Quantity of Interest: Cumulative human infections
                H_s_initial = u0[DengueModel.IX_Hₛ]
                H_s_final   = sol.u[end][DengueModel.IX_Hₛ]
                cumulative_infections = H_s_initial - H_s_final
                results[i] = max(0.0, cumulative_infections)
            else
                results[i] = 0.0
            end

        catch e
            # Error occurred - assign zero
            # (printing from many threads can get noisy; keep as-is but rare)
            if i % 100 == 0
                println("  Warning: Error in sample $i: $e")
            end
            results[i] = 0.0
        end
    end

    println("  ✓ Completed $num_samples simulations")
    return results
end

# 4. Run Sobol Analysis
println("\n" * "="^60)
println("EXECUTING SOBOL ANALYSIS")
println("="^60)

bounds = [(lower_bounds[i], upper_bounds[i]) for i in 1:length(lower_bounds)]

println("\nSampling scheme: Sobol sequence")
println("Number of samples: 500")

start_time = time()
res = gsa(sobol_interface, Sobol(), bounds; samples=50)
elapsed = time() - start_time

println("\n✓ Analysis completed in $(round(elapsed/60, digits=1)) minutes")

# 5. Extract Results
st_values = vec(collect(res.ST))
s1_values = vec(collect(res.S1))

# 6. Display Results
println("\n" * "="^60)
println("RESULTS: SENSITIVITY INDICES")
println("="^60)

println("\nFirst-Order Indices (S1) - Direct effect only:")
println("-"^60)
for i in 1:length(param_names)
    @printf("  %-8s: %.4f\n", param_names[i], s1_values[i])
end

println("\nTotal Sensitivity Indices (ST) - Total effect (with interactions):")
println("-"^60)
for i in 1:length(param_names)
    @printf("  %-8s: %.4f\n", param_names[i], st_values[i])
end

println("\nInteraction Effects (ST - S1):")
println("-"^60)
for i in 1:length(param_names)
    interaction = st_values[i] - s1_values[i]
    @printf("  %-8s: %.4f\n", param_names[i], interaction)
end

println("\nRanking by Total Sensitivity (ST):")
println("-"^60)
sorted_indices = sortperm(st_values, rev=true)
for (rank, idx) in enumerate(sorted_indices)
    @printf("  %d. %-8s (ST = %.4f)\n", rank, param_names[idx], st_values[idx])
end

# 8. Generate Plots
println("\n" * "="^60)
println("GENERATING PLOTS")
println("="^60)

p1 = bar(string.(param_names), st_values,
    title="Total Sensitivity Index (ST)",
    ylabel="Sensitivity",
    xlabel="Parameter",
    legend=false,
    color=:teal,
    xrotation=45,
    ylim=(0, maximum(st_values)*1.2),
    size=(800, 600)
)

p2 = bar(string.(param_names), s1_values,
    title="First-Order Sensitivity Index (S1)",
    ylabel="Sensitivity",
    xlabel="Parameter",
    legend=false,
    color=:coral,
    xrotation=45,
    ylim=(0, maximum(st_values)*1.2),
    size=(800, 600)
)

p_combined = plot(p1, p2, layout=(2,1), size=(800, 1000),
    plot_title="Dengue Model Sensitivity Analysis"
)

mkpath("outputs")
savefig(p_combined, "outputs/sobol/dengue_sobol_analysis.png")
println("  ✓ Plot saved: outputs/dengue_sobol_analysis.png")

# 9. Save Numerical Results
using DataFrames, CSV

results_df = DataFrame(
    Parameter = param_names,
    S1 = s1_values,
    ST = st_values,
    Interaction = st_values .- s1_values,
    Rank = [findfirst(==(i), sorted_indices) for i in 1:length(param_names)]
)

CSV.write("outputs/dengue_sobol_results.csv", results_df)
println("  ✓ Results saved: outputs/dengue_sobol_results.csv")

println("\n" * "="^60)
println("ANALYSIS COMPLETE")
println("="^60)
