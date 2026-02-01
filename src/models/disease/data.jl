using CSV
using DataFrames
using Dates
using Statistics

export get_dengue_data

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
function get_dengue_data(; filename::AbstractString = "../data/dengue_cases-2010_2022.csv",
                          mean::Bool = true)::DataFrame
    path = _relpath(filename)
    df = CSV.read(path, DataFrame)

    @assert "dt_sin_pri" in names(df) "Expected a 'dt_sin_pri' column in $path"
    df.dt_sin_pri = Date.(df.dt_sin_pri)
    sort!(df, :dt_sin_pri)

    if mean
        df = rolling_mean_df(df, :dt_sin_pri, 7)
    end

    return df
end