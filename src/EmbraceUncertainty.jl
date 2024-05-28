module EmbraceUncertainty

using Arrow
using CSV
using DataFrames
using Dates
using Downloads
using Markdown
using MixedModels
using PooledArrays
using Scratch
using ZipFile

const CACHE = Ref("")
const MMDS = String[]
const ML_LATEST_URL = "https://files.grouplens.org/datasets/movielens/ml-latest.zip"

_file(x) = joinpath(CACHE[], string(x, ".arrow"))

function __init__()
    CACHE[] = @get_scratch!("data")
    append!(MMDS, MixedModels.datasets())
end

clear_scratchspaces!() = Scratch.clear_scratchspaces!(@__MODULE__)

function extract_csv(zipfile, fname; kwargs...)
    file = only(filter(f -> endswith(f.name, fname), zipfile.files))
    return CSV.read(file, DataFrame; delim=',', header=1, kwargs...)
end

const metadata = Dict{String,String}("url" => ML_LATEST_URL)

function create_arrow(fname, df)
    arrowfile = _file(splitext(basename(fname))[1])
    Arrow.write(arrowfile, df; compress=:lz4, metadata)
    return arrowfile
end

const GENRES = ["Action", "Adventure", "Animation",
                "Children", "Comedy", "Crime",
                "Documentary", "Drama",
                "Fantasy", "Film-Noir",
                "Horror",
                "IMAX",
                "Musical", "Mystery",
                "Romance",
                "Sci-Fi",
                "Thriller",
                "War", "Western"]

function load_quiver()
    @info "Downloading data"
    quiver = String[]
    open(Downloads.download(ML_LATEST_URL), "r") do io
        zipfile = ZipFile.Reader(io)
        @info "Extracting and saving ratings"
        ratings = extract_csv(zipfile, "ratings.csv";
            types=[Int32, Int32, Float32, Int32],
            pool=[false, false, true, false],
        )
        push!(quiver, create_arrow("ratings.csv", ratings))
        @info "Extracting movies that are in the ratings table"
        movies = leftjoin!(
            leftjoin!(
                sort!(combine(groupby(ratings, :movieId), nrow => :nrtngs), :nrtngs),
                extract_csv(zipfile, "movies.csv"; types=[Int32,String,String], pool=false);
                on=:movieId,
            ),
            extract_csv(zipfile, "links.csv"; types=[Int32,Int32,Int32]);
            on=:movieId,
        )
        disallowmissing!(movies; error=false)
        movies.nrtngs = Int32.(movies.nrtngs)
        for g in GENRES
            setproperty!(movies, replace(g, "-" => ""), contains.(movies.genres, g))
        end
        select!(movies, Not("genres"))  # now drop the original genres column
        push!(quiver, create_arrow("movies.csv", movies))
        @info "Extracting and saving README"
        readme = only(filter(f -> endswith(f.name, "README.txt"), zipfile.files))
        open(joinpath(CACHE[], "README.txt"), "w") do io
            write(io, read(readme))
        end

        return nothing
    end

    return quiver
end

const OSF_IO_URIs = Dict{String,String}(
    "box" => "tkxnh",
    "elstongrizzle" => "5vrbw",
    "oxboys" => "cz6g3",
    "sizespeed" => "kazgm",
    "ELP_ldt_item" => "c6gxd",
    "ELP_ldt_subj" => "rqenu",
    "ELP_ldt_trial" => "3evhy",
    "movies" => "kvdch",
    "ratings" => "v73ym",
)

function osf_io_dataset(name::AbstractString)
    if haskey(OSF_IO_URIs, name)
        Downloads.download(
            string("https://osf.io/", OSF_IO_URIs[name], "/download"),
            _file(name),
        )
        return true
    end
    return false
end        

"""
    dataset(name::Union(Symbol, AbstractString))

Return as an `Arrow.Table` the dataset named `name`.
"""
dataset(name::Symbol) = dataset(string(name))
function dataset(name::AbstractString)
    name in MMDS && return MixedModels.dataset(name)
    f = _file(name)
    isfile(f) || osf_io_dataset(name) ||
        throw(ArgumentError("$(name) is not a dataset "))
    return Arrow.Table(f)
end

function readme()
    return Markdown.parse_file(joinpath(CACHE[], "README.txt"))
end

"""
    tagpad(v::AbstractVector{<:Integer}, tag::String="S"; pool::Bool=true)

Convert `v` to a vector of strings prepended with `tag` and padded to a constant string length.
If `pool` is `true`, the default, the resulting vector of strings is converted to a `PooledArray`.

The reason for padding the numeric strings is so that the strings sort lexicographically in the
same order as the original numeric values.

The single-argument version, e.g. `tagpad(:I)`, returns a partially-applied function that can be used
in a `transform` or `select` call.

```@example
show(tagpad(repeat(1:10, inner=2)))
```
"""
function tagpad(v::AbstractVector{<:Integer}, tag::String="S"; pool::Bool=true)
    tagged = string.(tag, lpad.(v, maximum(ndigits, v), '0'))
    return pool ? PooledArray(tagged; signed=true, compress=true) : tagged
end

tagpad(v::AbstractVector{<:Integer}, tag; pool::Bool=true) = tagpad(v, string(tag); pool)

tagpad(tag) = Base.Fix2(tagpad, string(tag))

"""
    age_at_event(edate::Dates.TimeType, dob::Dates.TimeType)

Return the age in years at `edate` for a person born on `dob`.
"""
function age_at_event(edate::TimeType, dob::TimeType)
    (ey, em, ed) = yearmonthday(edate)
    (by, bm, bd) = yearmonthday(dob)
    return (ey - by) - (em < bm | (em == bm & ed < bd))
end

export GENRES,
    age_at_event,
    tagpad


end # module EmbraceUncertainty
