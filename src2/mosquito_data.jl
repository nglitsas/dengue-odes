# src2/mosquito_data.jl
# Data loading + preprocessing utilities for the MosquitoCapture sub-model.
#
# This is a Julia translation of the Python preprocessing helpers:
# - trap data loader (monthly aggregation, odd months only, trapped = dead+alive females)
# - dengue data loader (optional 7-day moving average)
# - weather data loader (optional date slicing)
# - temperature extraction helper
#
# IMPORTANT INTEGRATION NOTE:
# - Paths are resolved relative to this file (`@__DIR__`) so your module works regardless
#   of where you run Julia from.
# - Filenames are kept identical to the Python defaults.

using CSV
using DataFrames
using Dates
using Statistics

# -----------------------------
# Path helper (robust integration)
# -----------------------------

"""
    _relpath(rel)

Resolve `rel` relative to this file's directory. This prevents failures when the
current working directory differs from the project layout.
"""
_relpath(rel::AbstractString) = normpath(joinpath(@__DIR__, rel))

# -----------------------------
# Core loaders (mirroring Python)
# -----------------------------

"""
    get_trap_data(; filename="../data/mosq_aaeg_trap-2017_2022.csv", replace=true)

Loads mosquito trap capture data.

Pipeline (mirrors the Python function):
1) Reads CSV and parses the `date` column.
2) Aggregates to monthly totals (sum).
3) Keeps only odd months (Jan, Mar, May, Jul, Sep, Nov).
4) Creates `trapped = m_aaeg_f_m + m_aaeg_f_v`.
5) Sets each timestamp to the first day of the month.
6) If `replace=true`, replaces the first row where `trapped == 0` with the
   integer-average of the neighboring rows (previous and next).

Returns: DataFrame with at least `date` and `trapped` plus aggregated columns.
"""
function get_trap_data(; filename::AbstractString = "../data/mosq_aaeg_trap-2017_2022.csv",
                        replace::Bool = true)::DataFrame
    path = _relpath(filename)
    df = CSV.read(path, DataFrame)

    @assert "date" in names(df) "Expected a 'date' column in $path"
    df.date = Date.(df.date)
    sort!(df, :date)

    # Group by year-month to emulate pandas resample('1M').sum()
    df.year  = year.(df.date)
    df.month = month.(df.date)

    g = groupby(df, [:year, :month])

    # Sum numeric columns (skip year/month themselves)
    numeric_cols = Symbol[]
    for c in names(df)
        cs = Symbol(c)
        eltype(df[!, cs]) <: Number || continue
        cs in (:year, :month) && continue
        push!(numeric_cols, cs)
    end

    agg = combine(g, numeric_cols .=> (x -> sum(skipmissing(x))) .=> numeric_cols)

    # Month timestamps -> first day of month (Python maps month-end to day=1)
    agg.date = Date.(agg.year, agg.month, 1)

    # Keep odd months only (drop even months)
    even_months = Set([2, 4, 6, 8, 10, 12])
    agg = agg[.!in.(agg.month, Ref(even_months)), :]

    # trapped = m_aaeg_f_m + m_aaeg_f_v
    @assert ("m_aaeg_f_m" in names(agg)) && ("m_aaeg_f_v" in names(agg)) """
    Expected columns 'm_aaeg_f_m' and 'm_aaeg_f_v' in aggregated trap data.
    """
    agg.trapped = agg.m_aaeg_f_m .+ agg.m_aaeg_f_v

    sort!(agg, :date)

    # Optional replacement for first trapped==0 (same behavior as Python)
    if replace
        idxs = findall(==(0), agg.trapped)
        if !isempty(idxs)
            idx = first(idxs)
            if 1 < idx < nrow(agg)
                for c in names(agg)
                    cs = Symbol(c)
                    eltype(agg[!, cs]) <: Number || continue
                    agg[idx, cs] = Int(round((agg[idx - 1, cs] + agg[idx + 1, cs]) / 2))
                end
                # recompute trapped (since numeric cols were overwritten)
                agg.trapped = agg.m_aaeg_f_m .+ agg.m_aaeg_f_v
            else
                @warn "Found trapped==0 on a boundary row; cannot replace using neighbors."
            end
        else
            println("No zero values found in trapped")
        end
    end

    return agg
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


"""
    get_mos_trap(; filename="../../../data/mosq_trapped_model.csv")

Loads trap data already cleaned and ready to be used in ODE fitting.

This file is assumed to have been exported from pandas with an index column named
"Unnamed: 0".

Returns: DataFrame with a `date::Date` column.
"""
function get_mos_trap(; filename::AbstractString = "../data/mosq_trapped_model.csv")::DataFrame
    path = _relpath(filename)
    df = CSV.read(path, DataFrame)

    @assert "Unnamed: 0" in names(df) "Expected an 'Unnamed: 0' column in $path"
    df.date = Date.(df."Unnamed: 0")
    select!(df, Not("Unnamed: 0"))
    sort!(df, :date)

    return df
end


"""
    parse_date(date_str)

Transforms date strings by adding leading zeros to single-digit day/month parts
to ensure date parsing is unambiguous.

Example:
- "1/2/2017" -> "01/02/2017"
"""
function parse_date(date_str::AbstractString)::String
    parts = split(date_str, '/')
    padded = [length(p) == 1 ? "0" * p : p for p in parts]
    return join(padded, '/')
end


"""
    get_weather_data(; filename="../weather-2010_2022.csv", ini_date=nothing, end_date=nothing)

Loads climate data that has already been preprocessed and saved.

Assumes the CSV contains a date column named "Data".

If `ini_date` and `end_date` are provided as strings (YYYY-mm-dd or parseable),
the data is filtered to that inclusive interval.

Returns: DataFrame with `date::Date`.
"""
function get_weather_data(; filename::AbstractString = "../data/weather-2010_2022.csv",
                           ini_date = nothing,
                           end_date = nothing)::DataFrame
    path = _relpath(filename)
    df = CSV.read(path, DataFrame)

    @assert "Data" in names(df) "Expected a 'Data' column in $path"
    df.date = Date.(df.Data)
    select!(df, Not("Data"))
    sort!(df, :date)

    if (ini_date isa AbstractString) && (end_date isa AbstractString)
        d0 = Date(ini_date)
        d1 = Date(end_date)
        df = df[(df.date .>= d0) .& (df.date .<= d1), :]
    end

    return df
end


"""
    get_temp(start_date, end_date)

Returns a vector with the mean temperature in a given time interval.

Parameters:
- start_date: string in "YYYY-mm-dd" format (or parseable by Date)
- end_date: string in "YYYY-mm-dd" format (or parseable by Date)

Returns:
- Vector of values from the `temp_mean-celsius` column.
"""
function get_temp(start_date::AbstractString, end_date::AbstractString)
    df = get_weather_data()
    @assert "temp_mean-celsius" in names(df) "Expected column 'temp_mean-celsius' in weather data"

    d0 = Date(start_date)
    d1 = Date(end_date)

    sub = df[(df.date .>= d0) .& (df.date .<= d1), :]
    return collect(sub[!, "temp_mean-celsius"])
end


# -----------------------------
# Integration-facing entry point
# -----------------------------

"""
    load_data(; source=:model, replace=true, mean=true)

Primary entry point used by the MosquitoCapture module.

- `source=:model`   -> loads "../../../data/mosq_trapped_model.csv" (already cleaned for fitting)
- `source=:raw`     -> loads "../mosq_aaeg_trap-2017_2022.csv" and applies the monthly/odd-month pipeline

Also provides access to dengue and weather data (returned as a NamedTuple).

Returns a NamedTuple:
    (trap = <DataFrame>, dengue = <DataFrame>, weather = <DataFrame>)

You can ignore dengue/weather if the mosquito-only fitting doesn’t need them yet.
"""
function load_data(; source::Symbol = :model,
                    replace::Bool = true,
                    mean::Bool = true)

    trap =
        source === :model ? get_mos_trap() :
        source === :raw   ? get_trap_data(replace=replace) :
        error("Unknown source=$source. Use :model or :raw")

    dengue = get_dengue_data(mean=mean)
    weather = get_weather_data()

    return (trap = trap, dengue = dengue, weather = weather)
end


# -----------------------------
# Internal helper: rolling mean
# -----------------------------

"""
    rolling_mean_df(df, date_col, window)

Applies rolling mean to numeric columns excluding `date_col`, and drops the first
`window-1` rows (pandas dropna behavior after rolling mean).
"""
function rolling_mean_df(df::DataFrame, date_col::Symbol, window::Int)::DataFrame
    @assert window >= 1 "window must be >= 1"
    n = nrow(df)
    n < window && return DataFrame()

    # numeric columns to smooth
    smooth_cols = Symbol[]
    for c in names(df)
        cs = Symbol(c)
        cs == date_col && continue
        eltype(df[!, cs]) <: Number || continue
        push!(smooth_cols, cs)
    end

    out = df[window:end, :]
    out[!, date_col] = df[window:end, date_col]

    for c in smooth_cols
        x = Float64.(df[!, c])
        y = Vector{Float64}(undef, n - window + 1)
        s = cumsum(x)
        y[1] = s[window] / window
        for i in (window + 1):n
            y[i - window + 1] = (s[i] - s[i - window]) / window
        end
        out[!, c] = y
    end

    return out
end
