module Data

using CSV
using DataFrames
using Dates
using Statistics

export get_trap_data, save_trap_data

# -----------------------------
# Path Helper
# -----------------------------
# Fixes paths relative to src/models/mosquito/data.jl
# We need to go up 3 levels: mosquito -> models -> src -> root -> data
const DATA_DIR = normpath(joinpath(@__DIR__, "../../../data"))

function _get_path(filename::String)
    return joinpath(DATA_DIR, filename)
end

# -----------------------------
# Core Loaders
# -----------------------------

"""
    get_trap_data(; filename="raw/mosq_aaeg_trap-2017_2022.csv", replace=true)

Loads raw trap data, aggregates to odd months, cleans it, and computes `mfai_obvs`.
"""
function get_trap_data(; filename::String = "raw/mosq_aaeg_trap-2017_2022.csv",
                        replace::Bool = true)::DataFrame
    path = _get_path(filename)
    df = CSV.read(path, DataFrame)

    @assert "date" in names(df) "Expected a 'date' column in $path"
    df.date = Date.(df.date)
    sort!(df, :date)

    # Group by year-month
    df.year  = year.(df.date)
    df.month = month.(df.date)
    g = groupby(df, [:year, :month])

    # Sum numeric columns
    numeric_cols = Symbol[]
    for c in names(df)
        cs = Symbol(c)
        eltype(df[!, cs]) <: Number || continue
        cs in (:year, :month) && continue
        push!(numeric_cols, cs)
    end

    agg = combine(g, numeric_cols .=> (x -> sum(skipmissing(x))) .=> numeric_cols)
    agg.date = Date.(agg.year, agg.month, 1)

    # Keep odd months only
    even_months = Set([2, 4, 6, 8, 10, 12])
    agg = agg[.!in.(agg.month, Ref(even_months)), :]

    # Calculate Total Trapped (Females)
    agg.trapped = agg.m_aaeg_f_m .+ agg.m_aaeg_f_v
    sort!(agg, :date)

    # Replace initial zeros
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
                agg.trapped = agg.m_aaeg_f_m .+ agg.m_aaeg_f_v
            end
        end
    end

    # --- Compute MFAI Observed (variable traps sampled) ---
    @assert :n_traps in Symbol.(names(agg)) "Expected column :n_traps to compute mfai_obvs"
    agg.mfai_obvs = Float64.(agg.trapped) ./ Float64.(agg.n_traps)

    return agg
end


"""
    save_trap_data(; in_file="raw/mosq_aaeg_trap-2017_2022.csv",
                    out_file="processed/mosq_trapped_model.csv",
                    replace=true)

Runs `get_trap_data` and writes the processed table to `data/processed/...`.
Returns the processed DataFrame.
"""
function save_trap_data(; in_file::String = "raw/mosq_aaeg_trap-2017_2022.csv",
                         out_file::String = "processed/mosq_trapped_model.csv",
                         replace::Bool = true)::DataFrame
    df = get_trap_data(; filename=in_file, replace=replace)

    out_path = _get_path(out_file)
    mkpath(dirname(out_path))
    CSV.write(out_path, df)

    return df
end

end # module
