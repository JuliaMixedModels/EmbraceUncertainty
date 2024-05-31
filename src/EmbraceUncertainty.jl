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
using SHA
using TypedTables
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

const OSF_IO_URIs = 
    CSV.read(
        IOBuffer(
"""
dsname,filename,version,sha2
box,tkxnh,1,4ac2038735c5286ca0d6a706b4feef5b34bd93560bc4cabc7223addf2366e4c0
elstongrizzle,5vrbw,1,f33e08ad5a91ab5dd9889953dbb97cc52a299fdac9db8b37ec4f87ad2dacadbd
oxboys,cz6g3,1,079f46e404e43c1848a65d3a4e5b6a14b2cb17565c77819d7ee3effe72f5ebd0
sizespeed,kazgm,1,d8eddc7f26928def4bff0e2c3637a90ff45fe589e3ffd5ab4818c8c1f52cf87e
ELP_ldt_item,c6gxd,1,f851910e1435659ca662ad49cfb8deb6b7cf287e4ce4969103dba11b32ab2e6c
ELP_ldt_subj,rgenu,2,d9c88915681b64fc9db975f9bb2d6f402058fee5cb35887f9de7d07776efdd56
ELP_ldt_trial,3evhy,2,57a83679f8f458b1d9bb56e05099a556ffdaf15af67464e9ee608c844fc4fa9c
movies,kvdch,1,c8fa488be74c368530f38de3c1d3511d2182e3fa07f357d2fa09121adc1cc964
ratings,v73ym,1,1d6466b8fd8da2942881c83941077cd2b6c9f7a03fa2d71072a5cd7aaa4ef560
fggk21,vwecy,1,0fa959f095f8b92135496b6f8c8a8b5a3e896e8875f0ba6928bd074559d8a796
fggk21_Child,c2fmn,1,61c91e00336e6f804e9f6b86986ebb4a14561cc4908b3a21cb27c113d2b51a5c
fggk21_Score,7fqx3,1,99d73ee705aaf5f4ee696eadbba992d0113ba6f467ce337a62a63853e4617400
kkl15,p8cea,2,90d7bb137c8613d7a15c8597c461aee7c7cb0f0989a07c80fc93e1fbe2e5c156
kwdyz11,4cv52,3,2fa23aa8aa25e1adb10183c8d29646ae0d19d6baef9d711c9906f7fa1b225571
"""
        ),
        Table,
    )

function osf_io_dataset(name::AbstractString)
    r = only(filter(==(name) ∘ getproperty(:dsname), OSF_IO_URIs))
    fn = _file(name)
    Downloads.download(
        string("https://osf.io/", r.filename, "/download?version=", r.version),
        fn,
    )
    if bytes2hex(open(sha2_256, fn)) ≠ r.sha2
        @error("sha2 of downloaded file doesn't match entry in table")
    end
    return true
end        

"""
    dataset(name::Union(Symbol, AbstractString); reload::Bool=true)

Return as an `Arrow.Table` the dataset named `name`.

If `reload` is `true` the dataset will first be downloaded from the osf.io site, even if a current copy exists.
"""
dataset(name::Symbol; reload::Bool=false) = dataset(string(name); reload)
function dataset(name::AbstractString; reload::Bool=false)
    name in MMDS && return MixedModels.dataset(name)
    f = _file(name)
    if reload | !isfile(f)
        if !osf_io_dataset(name)
            throw(ArgumentError("$(name) is not a dataset "))
        end
    end
    return Arrow.Table(f)
end

datasets() = OSF_IO_URIs

function readme()
    return Markdown.parse_file(joinpath(CACHE[], "README.txt"))
end

"""
    tagpad(v::AbstractVector{<:Integer}, ndig::Integer, tag::String="S"; pool::Bool=true)
    tagpad(v::AbstractVector{<:Integer}, tag::String="S"; pool::Bool=true)

Convert `v` to a vector of strings prepended with `tag` and padded to a constant string length of `ndig`,
which is evaluated as `maximum(ndigits, v)`, if not provided.

If `pool` is `true`, the default, the resulting vector of strings is converted to a `PooledArray`.

The reason for padding the numeric strings is so that the strings sort lexicographically in the
same order as the original numeric values.

The single-argument version, e.g. `tagpad(:I)`, returns a partially-applied function that can be used
in a `transform` or `select` call.

```@example
show(tagpad(repeat(1:10, inner=2)))
```
"""
function tagpad(v::AbstractVector{<:Integer}, ndig::Integer, tag::AbstractString="S"; pool::Bool=true)
    tagged = string.(tag, lpad.(v, ndig, '0'))
    return pool ? PooledArray(tagged; signed=true, compress=true) : tagged
end

function tagpad(v::AbstractVector{<:Integer}, tag::AbstractString="S"; pool::Bool=true)
    return tagpad(v, maximum(ndigits, v), tag; pool)
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
