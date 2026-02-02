module Constants

export N_TRAPS, N_HOUSEHOLDS, ALPHA, K, CAPTURE_RATE, POPULATION

# Foz do Iguaçu parameters (Table 1)
const N_TRAPS = 2412.0
const N_HOUSEHOLDS = 102751.0
const ALPHA = 0.0007
const K = 0.5                # Fraction of females
const POPULATION = 256088.0
# Derived constant for the ODE
# Rate at which mosquitoes are removed from M and added to Trapped
const CAPTURE_RATE = ALPHA * (N_TRAPS / N_HOUSEHOLDS)

end