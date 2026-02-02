module Data
using CSV
using DataFrames
using Dates
using Statistics
using RollingFunctions # Might need to add this package: ] add RollingFunctions

export get_dengue_data

# --- ADDED HELPER: Rolling Mean ---
function rolling_mean_df(df::DataFrame, date_col::Symbol, window::Int)
    # create a copy to avoid modifying original
    res = copy(df)
    
    # identify numeric columns to smooth
    cols_to_smooth = names(res, Number)
    
    # Apply rolling mean to each numeric column
    for c in cols_to_smooth
        # use runmean from RollingFunctions or manual calculation
        res[!, c] = runmean(res[!, c], window)
    end
    
    # Filter out the first few rows where the rolling mean is incomplete
    # (The first 'window-1' rows usually have artifacts)
    return res[window:end, :]
end

"""
    get_dengue_data(; filename="../dengue_cases-2010_2022.csv", mean=true)

Loads dengue data previously cleaned.

Returns time series of reported/probable/lab-confirmed cases (depending on columns
in your CSV).

Parameters:
- `mean`: If true, applies a 7-day moving average and drops the initial window,
          mirroring pandas rolling(window=7).mean().dropna().

Returns: DataFrame with `dt_sin_pri::Date` and remaining columns.
"""
function get_dengue_data(; filename::AbstractString = "data/raw/dengue_cases-2010_2022.csv",
                           mean::Bool = true)::DataFrame
    path = abspath(filename)
    df = CSV.read(path, DataFrame)

    @assert "dt_sin_pri" in names(df) "Expected a 'dt_sin_pri' column in $path"
    df.dt_sin_pri = Date.(df.dt_sin_pri)
    sort!(df, :dt_sin_pri)

    if mean
        df = rolling_mean_df(df, :dt_sin_pri, 7)
    end

    return df
end

end # module