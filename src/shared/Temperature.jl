module Temperature

using DataFrames
using CSV
using Dates
using DataInterpolations
using Statistics

using ..TimeUtil

export get_weather_data, get_temperature_interpolator

const DATA_DIR       = normpath(joinpath(@__DIR__, "../../data"))
const RAW_PATH       = joinpath(DATA_DIR, "raw/weather-2010_2022.csv")
const PROCESSED_PATH = joinpath(DATA_DIR, "processed/temperature_2010_2022.csv")

function get_weather_data(; force_process::Bool=false)::DataFrame
    mkpath(dirname(PROCESSED_PATH))

    if isfile(PROCESSED_PATH) && !force_process
        println("Loading processed weather data from: $PROCESSED_PATH")
        df = CSV.read(PROCESSED_PATH, DataFrame)
        df.date = Date.(df.date)
        return df  # ✓ Now includes t_vals column from CSV
    end

    println("Processing raw weather data from: $RAW_PATH")
    if !isfile(RAW_PATH)
        error("Raw weather file not found at: $RAW_PATH")
    end

    # 1. Read Raw Data
    df = CSV.read(RAW_PATH, DataFrame)
    
    # 2. Rename Columns
    rename!(df, "Data" => "date", "temp_mean_celsius" => "temp")
    df.date = Date.(df.date)
    
    # 3. Sort
    sort!(df, :date)
    
    # 4. Add t_vals column (days since 2000)
    df.t_vals = Float64.(date_to_t.(df.date))
    
    select!(df, :date, :t_vals, :temp)  # ✓ Include t_vals

    # 5. Save Processed Data
    CSV.write(PROCESSED_PATH, df)
    println("Saved processed weather data to: $PROCESSED_PATH")

    return df
end

function get_temperature_interpolator(weather_df::DataFrame)
    # 1. Sort by t_vals
    sort!(weather_df, :t_vals)

    # 2. Filter missing
    clean_df = dropmissing(weather_df, [:t_vals, :temp])

    if nrow(clean_df) < nrow(weather_df)
        println("Warning: Dropped $(nrow(weather_df) - nrow(clean_df)) rows with missing data.")
    end

    # 3. Extract vectors
    clean_df.t_vals = Float64.(clean_df.t_vals)
    clean_df.temp = Float64.(clean_df.temp)
    sort!(clean_df, :t_vals)
    
    @assert !any(isnan, clean_df.temp) "Weather data contains NaNs!"
    @assert !any(isnan, clean_df.t_vals)    "Time values contain NaNs!"

    # 4. Create Interpolator
    return LinearInterpolation(clean_df.temp, clean_df.t_vals)
end
end 