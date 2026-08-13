# Regenerate the `sizespeed` table and the `optsums/mvm*.json` files for
# @sec-largescaleobserved.
#
# Run this once, from the repository root, after changing the MovieLens data or
# the version of MixedModels.jl used to fit these models:
#
#     julia --startup-file=no --project -t auto scripts/movielens_sizespeed.jl
#
# It writes `src/sizespeed.csv` and one `optsums/mvmXXuYY.json` per cutoff pair.
# Both are tracked in git, because refitting is far too expensive to do while
# rendering the book: `restoreoptsum!` re-evaluates the objective and rejects a
# file that does not belong to the data at hand, so these outputs must be
# regenerated together with any change to the ratings table.
#
# Resource requirements: the `mc = 1` fits carry an 84,432 x 84,432 lower
# triangle in the [2,2] block of `L`, which is about 27 GiB even in the
# rectangular full packed format that MixedModels.jl v5.8.0 and later use by
# default.  Budget well over 32 GiB of RAM and expect the full sweep to take
# many hours.
#
# `nv` and `fittime` are not reproducible across machines --- they shift with
# the BLAS implementation and with trivial differences in floating point
# results --- but `evtime` is reasonably stable, as discussed in the chapter.

using CSV
using DataFrames
using EmbraceUncertainty: dataset, tagpad
using MixedModels
using ProgressMeter
using TypedTables

const MCUTOFFS = (1, 2, 5, 10, 15, 20, 50)
const UCUTOFFS = (20, 40, 80)

const FORM = @formula(rating ~ 1 + (1 | userId) + (1 | movieId))

optsumdir(paths::AbstractString...) =
  joinpath(@__DIR__, "..", "optsums", paths...)

"""
    ratingstable()

Return the MovieLens 32M ratings as a `Table` with the `nrtngs` and `urtngs`
columns that the cutoffs are applied to.

This mirrors the data preparation in @sec-largescaleobserved, except that the
movie titles are not needed here.
"""
function ratingstable()
  ratings = dataset(:ml32_ratings)
  df = DataFrame(;
    userId=ratings.userId,
    movieId=ratings.movieId,
    rating=Float32.(ratings.rating),
  )
  transform!(groupby(df, :movieId), nrow => :nrtngs)
  transform!(groupby(df, :userId), nrow => :urtngs)
  df.userId = tagpad(df.userId, "U")
  df.movieId = tagpad(df.movieId, "M")
  return Table(df)
end

"""
    fitcutoffs(data, mcutoff, ucutoff)

Fit the model to `data` trimmed to at least `mcutoff` ratings per movie and
`ucutoff` ratings per user, save its `optsum`, and return a `NamedTuple` of the
row to record in `sizespeed`.
"""
function fitcutoffs(data, mcutoff::Integer, ucutoff::Integer)
  trimmed = filter(
    r -> (r.nrtngs ≥ mcutoff) & (r.urtngs ≥ ucutoff),
    data,
  )
  model = LinearMixedModel(FORM, trimmed)
  model.optsum.initial .= 0.5
  fittime = @elapsed fit!(model)

  saveoptsum(
    optsumdir(
      "mvm$(lpad(mcutoff, 2, '0'))u$(lpad(ucutoff, 2, '0')).json",
    ),
    model,
  )

  # time a single objective evaluation at the optimum
  θ = copy(model.optsum.final)
  objective(updateL!(setθ!(model, θ)))    # warm up
  evtime = @elapsed objective(updateL!(setθ!(model, θ)))

  # `reterms` are ordered by decreasing number of levels, so look the grouping
  # factors up by name rather than by position
  nlevels = Dict(
    string(f) => length(t.levels) for
    (f, t) in zip(fnames(model), model.reterms)
  )

  # the integer columns are written as plain digits; `dataset(:sizespeed)`
  # reads them back with `downcast=true`, which picks the narrowest type
  return (;
    mc=mcutoff,
    uc=ucutoff,
    nratings=length(trimmed),
    nusers=nlevels["userId"],
    nmvie=nlevels["movieId"],
    modelsz=Base.summarysize(model) / 2^30,
    # blocks are stored in a one-dimensional array: 1 = [1,1], 2 = [2,1],
    # 3 = [2,2], which is the block that dominates the memory footprint
    L22sz=Base.summarysize(model.L[3]) / 2^30,
    nv=model.optsum.feval,
    fittime,
    evtime,
  )
end

function main()
  @info "Building the ratings table"
  data = ratingstable()

  rows = NamedTuple[]
  @showprogress showspeed=true for ucutoff in UCUTOFFS, mcutoff in MCUTOFFS
    @info "Fitting mc = $mcutoff, uc = $ucutoff"
    row = fitcutoffs(data, mcutoff, ucutoff)
    @info "  " row...
    push!(rows, row)
    GC.gc()
  end

  target = joinpath(@__DIR__, "..", "src", "sizespeed.csv")
  CSV.write(target, DataFrame(rows))
  @info "Wrote $target"
  return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
