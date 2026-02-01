module Shared

# --- Placeholders (Commented out because files are empty) ---
# When you write code in these files later, uncomment these lines.

include("TimeUtil.jl") 
using .TimeUtil
export TimeUtil

# --- Active Modules (Files that have code) ---
include("Temperature.jl")
using .Temperature
export Temperature

# include("IO.jl")
# using .IO
# export IO

# include("Losses.jl")
# using .Losses
# export Losses

# include("Params.jl")
# using .Params
# export Params

# include("Solve.jl")
# using .Solve
# export Solve

# include("Forcing.jl")
# using .Forcing
# export Forcing

end # module