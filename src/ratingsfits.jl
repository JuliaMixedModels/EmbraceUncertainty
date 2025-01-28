# this script is used to generate the various model fits for the movielens data
# and cache the results, including timing information
using Arrow
using DataFrames
using EmbraceUncertainty: EmbraceUncertainty as EU
using MixedModels

using MixedModels: ProgressMeter

rootpath() = dirname(joinpath(@__DIR__))

datadir(paths::AbstractString...) = joinpath(rootpath(), "data", paths...)
optsumdir(paths::AbstractString...) = joinpath(rootpath(), "optsums", paths...)

mkpath(datadir())
mkpath(optsumdir())

# make sure we're using the 32M dataset
EU.load_quiver(EU.ML_32M_URL)
ratings = DataFrame(Arrow.Table(joinpath(EU.CACHE[], "ratings.arrow")))
disallowmissing!(
    leftjoin!(
        leftjoin!(
            select!(ratings, Not(:timestamp)),
            combine(groupby(ratings, :movieId), nrow => :nrtngs);
            on=:movieId,
        ),
        combine(groupby(ratings, :userId), nrow => :urtngs);
        on=:userId,
    ),
)

function dooptsum(
    mcutoff::Integer,
    ucutoff::Integer;
    data=ratings,
    form=@formula(rating ~ 1 + (1|userId) + (1|movieId)),
    fitlog=true,
    initial=[0.5, 0.5],
    initial_step=[0.1,0.1],
    progress=false,
)
    mvm = LinearMixedModel(form, filter(r -> r.nrtngs ≥ mcutoff && r.urtngs ≥ ucutoff, data))
    mvm.optsum.initial = initial
    mvm.optsum.initial_step = initial_step
    mvusr = size(mvm.L[2])
    println()
    fittm = @timed fit!(mvm; fitlog, progress)
    saveoptsum(optsumdir("mvm$(lpad(mcutoff, 2, '0'))u$(lpad(ucutoff, 2, '0')).json"), mvm)
    evtm = @timed objective(updateL!(setθ!(mvm, mvm.θ)))
    open(datadir("sizespeed.csv"), "a") do io
        join(
            io,
            (
                mcutoff,
                ucutoff,
                size(mvm.Xymat, 1),       # number of observations
                last(mvusr),              # number of users
                first(mvusr),             # number of movies
                Float32(Base.summarysize(mvm)/2^30),      # size of model object (GiB)
                Float32(Base.summarysize(mvm.L[3])/2^30), # size of L22 (GiB)
                length(mvm.optsum.fitlog),# total number of objective function evaluations
                Float32(fittm.time),
                Float32(evtm.time),
            ),
            ',',
        )
        println(io)
    end
    return mvm
end

function ratingsoptsum(
    mcutoff::Integer,
    ucutoff::Integer;
    data=ratings,
    form=@formula(rating ~ 1 + (1|userId) + (1|movieId)),
)
    optsumfnm = optsumdir("mvm$(lpad(mcutoff, 2, '0'))u$(lpad(ucutoff, 2, '0')).json")
    isfile(optsumfnm) || throw(ArgumentError("File $optsumfnm is not available"))
    return restoreoptsum!(
        LinearMixedModel(form, filter(r -> r.nrtngs ≥ mcutoff && r.urtngs ≥ ucutoff, data)),
        optsumfnm,
    )
end
usteps = [5, 10, 20, 40, 80]
msteps = [1, 2, 5, 10, 15, 20, 50]
p = ProgressMeter.Progress(length(usteps) * length(msteps); showspeed=true)
for u in usteps, m in msteps
    GC.gc()
    GC.gc()
    dooptsum(m, u)
    GC.gc()
    ProgressMeter.next!(p)
end
ProgressMeter.finish!(p)
