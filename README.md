# Dengue-ODES

Temperature-dependent ordinary differential equation (ODE) models for **mosquito population dynamics** and **dengue transmission**, inspired by the methodology in Rauh et al. (2025) for Foz do Iguaçu, Brazil.

The repository includes:
- A standalone mosquito capture model calibrated to trap data (MFAI)
- A coupled human–mosquito dengue transmission model
- Sequential calibration pipelines (mosquito → disease)
- Reproducible fitting, simulation, and reporting scripts

## Repository Structure
Dengue-ODES/
├── configs/               # Parameter bounds, fit windows, fixed params, paths
├── data/
│   ├── raw/               # Original trap, case, weather data
│   └── processed/         # Cleaned CSVs (mosq_trapped_model.csv, temperature_2010_2022.csv)
├── scripts/               # Entry-point executables
├── src/                   # Core Julia package (DengueODES)
│   ├── DengueODES.jl
│   ├── models/
│   │   ├── disease/       # Dengue full model
│   │   └── mosquito/      # Mosquito-only capture model
│   ├── observe/           # Observation models (MFAI, incidence)
│   └── shared/            # Utilities (temperature, time, forcing, etc.)
├── outputs/
│   ├── mosquito_fit/      # Fitted params, plots, stats
│   ├── disease_fit/
│   └── mosquito_sim/
├── test/
├── archive/               # Older / alternative implementations
├── Project.toml
├── Manifest.toml
└── README.md
text## Core Package Layout (src/DengueODES.jl)

All modeling logic is exposed via the `DengueODES` package.

### Shared Utilities (`src/shared/`)

- `TimeUtil.jl` — Date ↔ continuous time conversion
- `Temperature.jl` — Weather data → interpolation / forcing functions
- `Forcing.jl` — Temperature-dependent rates
- `Params.jl` — Parameter structs
- `IO.jl` — Data I/O helpers

### Observation Models (`src/observe/`)

- `MFAI.jl` — Theoretical mosquito-female adult index from trap model
- `Incidence.jl` — Reported dengue cases from human compartments

### Mosquito Model (`src/models/mosquito/`)

Temperature-sensitive aquatic + adult dynamics calibrated to trap MFAI.

| File              | Purpose                              |
|-------------------|--------------------------------------|
| `constants.jl`    | Fixed biological / site parameters   |
| `model_dynamics.jl` | ODE right-hand side (`CaptureModel!`) |
| `fitting.jl`      | Parameter estimation (C₀, bₖ, ϵ)     |
| `Report.jl`       | Diagnostics, plots, summary tables   |
| `Simulate.jl`     | Forward simulation                   |

### Disease Model (`src/models/disease/`)

Full SEIR-like human + mosquito transmission (not yet fully detailed in README).

## Multithreaded / Parallel Fitting (Mosquito Model)

The mosquito parameter fitting (`src/models/mosquito/fitting.jl`) uses **Latin Hypercube Sampling (LHS)** to explore parameter space efficiently.

Key features:
- **Parallel evaluation** of 800–2000+ samples via Julia's built-in multithreading (`@threads`)
- Each ODE solve runs independently → near-linear speedup on multi-core machines (e.g., 6–7× on 8 threads)
- Handles unstable solves gracefully (large penalty for blow-ups)
- Uses stiff-capable solvers (e.g., `Rodas5P`, `Rosenbrock23`) for better stability

To take full advantage:
```bash
JULIA_NUM_THREADS=8 julia scripts/run_mosquito_fit.jl
# or make it permanent: export JULIA_NUM_THREADS=8
```
# Scripts (scripts/)

Reproducible entry points.

| Script                  | Purpose                                              |
|-------------------------|------------------------------------------------------|
| `run_mosquito_fit.jl`   | Calibrate mosquito model (C₀, bₖ, ϵ) to trap data    |
| `run_mosquito_sim.jl`   | Forward simulation with fixed/fitted params          |
| `run_dengue_fit.jl`     | Calibrate transmission parameters (in progress)      |
| `run_all.jl`            | Full sequential pipeline                             |
| `run_pipeline.jl`       | Alternative / experimental workflow                  |

Example usage:
```bash
JULIA_NUM_THREADS=8 julia --project scripts/run_mosquito_fit.jl
```
## File Structure
### Data

data/raw/ — Original sources (do not modify)
data/processed/ — Cleaned inputs for fitting/simulation

### Outputs
Results saved under outputs/:

mosquito_fit/ — fitted_params.csv, fit_plot.png, fit_stats.csv, mfai_observed_vs_pred.csv, etc.
mosquito_sim/ — simulation trajectories and summaries

## Installation & Running

1. Clone the repo:
```Bash 
git clone <your-repo-url>
cd Dengue-ODES
```
2. Install dependencies:
```julia
julia --project
] instantiate
```
3. Run the mosquito fit (with parallelism):
```Bash
JULIA_NUM_THREADS=8 julia scripts/run_mosquito_fit.jl
```

## Calibration Approach
Mosquito stage — Nonlinear least squares via LHS + multithreading
Objective:
$[\argmin_{\theta} \sum_i \left( \text{MFAI}_i - \widehat{\text{MFAI}}_i(\theta) \right)^2]$
where $θ = [C₀, bₖ, ϵ]$ and $ϵ$ is a climate dependent carrying capacity for the mosquito population.
Disease stage — Sequential (uses mosquito equilibrium as input).

## Contact
For questions, collaboration, bug reports, or to discuss using/extending the code, please reach out:

Sophie Zhou: sophiezy@umich.edu
Nicholas Litsas: nglitsas@umich.edu

You can also open an issue on GitHub — we welcome feedback!


## Reference
Rauh CS, Araujo EC, Ganem F, Lana RM, Leandro AS, Martins CA, et al. (2025)
Assessing mosquito dynamics and dengue transmission in Foz do Iguaçu, Brazil through an enhanced temperature-dependent mathematical model.
PLOS ONE 20(9): e0330902.
https://doi.org/10.1371/journal.pone.0330902
text