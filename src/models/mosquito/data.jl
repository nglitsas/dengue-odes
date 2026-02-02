module Data

using CSV
using DataFrames
using Dates
using Statistics

# Import TimeUtil so trap data aligns with global weather time
using ...Shared.TimeUtil

export get_trap_data, save_trap_data

const DATA_DIR = normpath(joinpath(@__DIR__, "../../../data"))

function _get_path(filename::String)
    return joinpath(DATA_DIR, filename)
end

function get_trap_data(; filename::String = "raw/mosq_aaeg_trap-2017_2022.csv",
                         replace::Bool = true)::DataFrame
    path = _get_path(filename)
    df = CSV.read(path, DataFrame)

    @assert "date" in names(df) "Expected a 'date' column in $path"
    df.date = Date.(df.date)
    sort!(df, :date)

    # Aggregation
    df.year  = year.(df.date)
    df.month = month.(df.date)
    g = groupby(df, [:year, :month])

    numeric_cols = Symbol[]
    for c in names(df)
        cs = Symbol(c)
        eltype(df[!, cs]) <: Number || continue
        cs in (:year, :month) && continue
        push!(numeric_cols, cs)
    end

    agg = combine(g, numeric_cols .=> (x -> sum(skipmissing(x))) .=> numeric_cols)
    agg.date = Date.(agg.year, agg.month, 1)

    even_months = Set([2, 4, 6, 8, 10, 12])
    agg = agg[.!in.(agg.month, Ref(even_months)), :]

    if "m_aaeg_f_m" in names(agg) && "m_aaeg_f_v" in names(agg)
        agg.trapped = agg.m_aaeg_f_m .+ agg.m_aaeg_f_v
    else
        agg.trapped = agg.m_aaeg_f_m
    end
    sort!(agg, :date)

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
                if "m_aaeg_f_m" in names(agg) && "m_aaeg_f_v" in names(agg)
                    agg.trapped = agg.m_aaeg_f_m .+ agg.m_aaeg_f_v
                else
                    agg.trapped = agg.m_aaeg_f_m
                end
            end
        end
    end

    # --- CRITICAL CHANGE ---
    # We use TimeUtil.date_to_t. 
    # This means for 2017 data, t will be ~6200.0, NOT 0.0.
    # This matches the weather data timeline.
    agg.t = TimeUtil.date_to_t.(agg.date)

    if :n_traps in Symbol.(names(agg))
        agg.mfai_obvs = Float64.(agg.trapped) ./ Float64.(agg.n_traps)
    else
        agg.mfai_obvs = Float64.(agg.trapped) ./ 1.0
    end
    return agg
end

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