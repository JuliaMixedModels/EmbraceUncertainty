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
using LinearAlgebra
using MixedModels
using ProgressMeter
using TypedTables

const MCUTOFFS = (1, 2, 5, 10, 15, 20, 50)
const UCUTOFFS = (20, 40, 80)

const FORM = @formula(rating ~ 1 + (1 | userId) + (1 | movieId))

# Optimizer settings for the whole sweep.  These are deliberately looser than
# the MixedModels defaults (`ftol_rel = 1e-12`, `ftol_abs = 1e-8`,
# `optimizer = :LN_NEWUOA`): at this scale a single objective evaluation costs
# between 45 seconds and half an hour, so the last few digits of the
# convergence criterion are not worth the hours they take.
#
# They are also recorded in each saved `optsum` and are restored along with it,
# so every row of `sizespeed` must be generated with the same values --- mixing
# settings across the grid would make `nv` and `fittime` incomparable between
# rows.  The fallback fit in `ratingsoptsum` in @sec-largescaleobserved sets
# these to match.
const FTOL_REL = 1.0e-10
const FTOL_ABS = 1.0e-6
const OPTIMIZER = :LN_BOBYQA

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
  model.optsum.ftol_rel = FTOL_REL
  model.optsum.ftol_abs = FTOL_ABS
  model.optsum.optimizer = OPTIMIZER
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

"""
    provenance()

Return the lines describing how this sweep was produced, to be written as a
comment header on `src/sizespeed.csv`.

`modelsz` and `L22sz` are determined by how MixedModels.jl represents the
model --- the rectangular full packed format introduced in v5.8.0 halved the
`[2,2]` block --- while `fittime` and `evtime` depend on the BLAS and on how
many threads it was given.  None of that is recoverable from the numbers
themselves, and `Manifest.toml` is not tracked in this repository, so record it
alongside the table.
"""
function provenance()
  blas = basename(first(BLAS.get_config().loaded_libs).libname)
  return [
    "generated by scripts/movielens_sizespeed.jl on " *
    Libc.strftime("%Y-%m-%d", time()),
    "julia $VERSION, MixedModels $(pkgversion(MixedModels)), " *
    "NLopt $(pkgversion(MixedModels.NLopt))",
    "$blas with $(BLAS.get_num_threads()) BLAS threads, " *
    "$(Threads.nthreads()) julia threads",
    "$(Sys.CPU_NAME) x$(Sys.CPU_THREADS)",
    "ftol_rel $FTOL_REL, ftol_abs $FTOL_ABS, optimizer $OPTIMIZER",
  ]
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
  open(target, "w") do io
    for line in provenance()
      println(io, "# ", line)
    end
    # `append` keeps the comment lines above: a plain `CSV.write(io, ...)`
    # rewinds the stream and discards them
    CSV.write(io, DataFrame(rows); append=true, writeheader=true)
  end
  @info "Wrote $target"
  return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
