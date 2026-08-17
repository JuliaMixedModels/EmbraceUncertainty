using Arrow, CSV, DataFrames, MixedModels, MKL

function timefit(
    mc::Integer,
    uc::Integer,
    ratings::AbstractDataFrame;
    initial_step::Vector{<:AbstractFloat} = [0.5, 0.5],
    ftol_rel::AbstractFloat = 1.0e-10,
    backend::Symbol = :nlopt,
    optimizer::Symbol = :LN_BOBYQA,
    )
    mc = Int8(mc)
    uc = Int8(uc)
    df = ratings[(ratings.mnrtngs .≥ mc) .& (ratings.unrtngs .≥ uc), :]
    nratings = Int32(size(df, 1))
    nusers = Int32(length(unique(df.userId)))
    nmvie = Int32(length(unique(df.movieId)))
    model = LinearMixedModel(@formula(rating ~ 1 + (1|userId) + (1|movieId)), df)
    model.optsum.initial_step = initial_step
    model.optsum.ftol_rel = ftol_rel
    model.optsum.ftol_abs = 1.0e-6
    model.optsum.backend = backend
    model.optsum.optimizer = optimizer
    modelsz = Float32(Base.summarysize(model) / (2^30))  # size in GiB
    L22sz = Float32(Base.summarysize(model.L[3]) / (2^30))
    @info mc, uc, nratings, nusers, nmvie, modelsz, L22sz
    flush(stdout)
    fittime = Float32(@elapsed fit!(model; progress=isinteractive()))
    nv = Int8(length(model.optsum.fitlog))
    mnevtm = fittime / nv
    return (; mc, uc, nratings, nusers, nmvie, modelsz, L22sz, nv, fittime, mnevtm)
end

function mktbl()
    ratings = DataFrame(Arrow.Table("./data/ratings.arrow"))
    res =  @NamedTuple{
        mc::Int8,
        uc::Int8,
        nratings::Int32,
        nusers::Int32,
        nmvie::Int32,
        modelsz::Float32,
        L22sz::Float32,
        nv::Int8,
        fittime::Float32,
        mnevtm::Float32,
    }[]
    for mc in Int8.([1, 2, 5, 10, 15, 20, 50])
        for uc in Int8.([20, 40, 80])
            try
                sizespeedrow = timefit(mc, uc, ratings)
                @info sizespeedrow
                flush(stdout)
                push!(res, sizespeedrow)
            catch e
            finally
                CSV.write("./data/sizespeed.csv", res)
            end
        end
    end
end

(@main)(args) = (mktbl(); nothing)
