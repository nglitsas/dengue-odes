module TimeUtil

using Dates

export GLOBAL_START_DATE, date_to_t, t_to_date

# --- Global Reference (Epoch) ---
# We use Jan 1, 2010. This must be ON or BEFORE the earliest date in your weather data.
const GLOBAL_START_DATE = Date(2000, 1, 1)

"""
    date_to_t(date::Date) -> Float64

Converts a real-world Date to model time `t` (days since GLOBAL_START_DATE).
Example: 2010-01-02 -> 1.0
"""
function date_to_t(d::Date)
    return Float64(Dates.value(d - GLOBAL_START_DATE))
end

"""
    t_to_date(t::Float64) -> Date

Converts model time `t` back to a real-world Date.
"""
function t_to_date(t::Real)
    return GLOBAL_START_DATE + Day(round(Int, t))
end

end # module