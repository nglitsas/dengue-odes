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
        return df
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
    
    # 3. Sort & Select
    sort!(df, :date)
    select!(df, :date, :temp)

    # 4. Save Processed Data (No filtering/clamping)
    CSV.write(PROCESSED_PATH, df)
    println("Saved processed weather data to: $PROCESSED_PATH")

    return df
end

function get_temperature_interpolator(weather_df::DataFrame)
    sort!(weather_df, :date)
    t_vals = date_to_t.(weather_df.date)
    temp_vals = Float64.(weather_df.temp)
    return LinearInterpolation(temp_vals, t_vals)
end

end # module