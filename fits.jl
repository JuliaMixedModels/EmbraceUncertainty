using CategoricalArrays
using DataFrames
using EmbraceUncertainty: dataset
using MixedModels
using SparseArrays         # for the nnz function
using Statistics           # for the mean function
using TypedTables
optsumdir(paths::AbstractString...) =
  joinpath(@__DIR__, "optsums", paths...)
@isdefined(contrasts) || const contrasts = Dict{Symbol,Any}()
@isdefined(progress) || const progress = false
ratings = dataset(:ratings)
ratings =
  Table(getproperties(ratings, (:userId, :movieId, :rating)))
movies = Table(
  disallowmissing!(
    leftjoin!(
      combine(
        groupby(DataFrame(ratings), :movieId),
        :rating => mean => :mnrtng,
      ),
      DataFrame(dataset(:movies));
      on=:movieId,
    );
    error=false,
  ),
)
users = Table(
  combine(
    groupby(DataFrame(ratings), :userId),
    nrow => :urtngs,
    :rating => mean => :umnrtng,
  ),
)
ratings = Table(
  disallowmissing!(
    leftjoin!(
      leftjoin!(
        DataFrame(ratings),
        DataFrame(getproperties(movies, (:movieId, :nrtngs)));
        on=:movieId,
      ),
      DataFrame(getproperties(users, (:userId, :urtngs)));
      on=:userId,
    ),
  ),
)
function ratingsoptsum(
  mcutoff::Integer,
  ucutoff::Integer;
  data=ratings,
  form=@formula(rating ~ 1 + (1 | userId) + (1 | movieId)),
  contrasts=Dict(:movieId => Grouping(), :userId => Grouping()),
)
  optsumfnm = optsumdir(
    "mvm$(lpad(mcutoff, 2, '0'))u$(lpad(ucutoff, 2, '0')).json",
  )
  model = LinearMixedModel(
      form,
      filter(
        r -> (r.nrtngs ≥ mcutoff) & (r.urtngs ≥ ucutoff),
        data,
      );
      contrasts,
    )
  isfile(optsumfnm) && return restoreoptsum!(model, optsumfnm; rtol=1e-6)

  @warn "File $optsumfnm is not available, fitting model."
  model.optsum.initial .= 0.5
  fit!(model; thin=1)
  saveoptsum(optsumfnm, model)
  return model
end

for m in [5, 10, 15, 20, 50], u in [5, 10, 20]
    @info "" m u
    GC.gc()
    GC.gc()
    GC.gc()
    ratingsoptsum(m, u)
end
