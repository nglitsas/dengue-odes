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

    # 4. Save Processed Data
    CSV.write(PROCESSED_PATH, df)
    println("Saved processed weather data to: $PROCESSED_PATH")

    return df
end

function get_temperature_interpolator(weather_df::DataFrame)
    # 1. Sort
    sort!(weather_df, :date)

    # 2. Filter missing
    clean_df = dropmissing(weather_df, [:date, :temp])

    if nrow(clean_df) < nrow(weather_df)
        println("Warning: Dropped $(nrow(weather_df) - nrow(clean_df)) rows with missing temperature data.")
    end

    # 3. Convert time using TimeUtil (Global Time: Days since 2000)
    sort!(clean_df, :date)
    t_vals = Float64.(date_to_t.(clean_df.date))
    
    # 4. Extract Temps
    temp_vals = Float64.(clean_df.temp)
    
    @assert !any(isnan, temp_vals) "Weather data contains NaNs!"
    @assert !any(isnan, t_vals)    "Time values contain NaNs!"

    # 5. Create Interpolator (Fixed Typo)
    return LinearInterpolation(
        temp_vals, 
        t_vals, 
        extrapolation = ExtrapolationType.Constant
    )
end

end # module