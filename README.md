# Dengue-ODES

This repository implements temperature-dependent ordinary differential equation (ODE) models
for mosquito population dynamics and dengue transmission, following the methodology of
Rauh et al. (2025) for Foz do Iguaçu, Brazil.

The project contains:

- A mosquito-only capture model calibrated to trap data
- A full human–mosquito dengue transmission model
- A sequential calibration pipeline
- Reproducible simulation and fitting scripts


## Top Level Repository Structure

Dengue-ODES/
├── data/
├── configs/
├── scripts/
├── src/
├── outputs/
├── test/
├── LICENSE
├── Manifest.toml
├── Project.toml
└── README.md



## Source Code Layout

All modeling code is organized as a single Julia package.

src/
└── DengueODES.jl



### Core Utilities

src/Core/


Shared infrastructure used by both mosquito and disease models.

| File | Purpose |
|------|----------|
| Time.jl | Date and time indexing |
| Temperature.jl | Climate preprocessing |
| Forcing.jl | Temperature-dependent parameters |
| Params.jl | Parameter and result structs |
| Solve.jl | ODE solver wrappers |
| Losses.jl | Objective functions |
| IO.jl | Data loading and saving |


### Observation Models

src/Observe/


Maps model states to observed data.

| File | Purpose |
|------|----------|
| MFAI.jl | Trap counts to $MFAI_{theo}$ |
| Incidence.jl | Human states to reported cases |


### Mosquito Model

src/Models/Mosquito/


Implements the mosquito capture model.

| File | Purpose |
|------|----------|
| RHS.jl | ODE system |
| Init.jl | Initial conditions |
| Simulate.jl | Forward simulation |
| Fit.jl | Fit C0, bcap, epsilon |
| Report.jl | Diagnostics and plots |


### Disease Model

src/Models/Disease/


Implements the dengue transmission model.

| File | Purpose |
|------|----------|
| RHS.jl | ODE system |
| Init.jl | Initial conditions |
| Simulate.jl | Forward simulation |
| Fit.jl | Fit transmission parameters |
| Report.jl | Diagnostics and plots |


## Scripts

scripts/


Reproducible execution pipeline.

| Script | Purpose |
|--------|----------|
| run_mosquito_fit.jl | Calibrate mosquito model |
| run_disease_fit.jl | Calibrate disease model |
| run_all.jl | Full pipeline |


Example:

```bash
julia --project scripts/run_mosquito_fit.jl
Data
data/
├── raw/
└── processed/
Contains trap counts, dengue cases, and climate data.
Raw data should not be edited.

Outputs
outputs/
├── mosquito_fit/
└── disease_fit/
Stores fitted parameters, trajectories, and figures.

Modeling Pipeline
Stage 1: Mosquito Calibration
Uses trap data only.

Input:

Climate data

Trap counts

Fitted parameters:

C0

bcap

epsilon

Output:

MosquitoFitResult

Saved in outputs/mosquito_fit/

Stage 2: Disease Calibration
Uses mosquito outputs and case data.

Input:

Mosquito artifacts

Climate data

Dengue cases

Fitted parameters:

phi

ab, am, ah

Transmission parameters

Output:

DiseaseFitResult

Saved in outputs/disease_fit/

Parameter Management
All parameters are defined in:

src/Core/Params.jl
Fixed parameters come from literature.
Only sensitive parameters are estimated.

Calibration Method
Model fitting uses nonlinear least squares:

min_theta sum_i (Y_i - Yhat_i(theta))^2

Implemented using Levenberg–Marquardt or Optimization.jl.

Installation
Clone the repository:

git clone <repository-url>
cd Dengue-ODES
Instantiate environment:

julia --project
] instantiate
Running the Full Pipeline
julia --project scripts/run_all.jl
Testing
julia --project
] test
Extending the Codebase
Modify climate forcing:
src/Core/Temperature.jl

Add observation models:
src/Observe/

Change mosquito biology:
src/Models/Mosquito/RHS.jl

Add Bayesian inference:
Replace Fit.jl with Turing.jl workflows

Reference
Rauh et al. (2025)
Assessing mosquito dynamics and dengue transmission in Foz do Iguaçu
PLOS ONE 20(9): e0330902