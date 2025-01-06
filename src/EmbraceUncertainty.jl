module EmbraceUncertainty

using Arrow
using CSV
using DataFrames
using Dates
using Downloads
using Markdown
using MixedModelsDatasets
using PooledArrays
using Scratch
using SHA
using TypedTables
using ZipFile

const CACHE = Ref("")
const MMDS = String[]

function __init__()
    CACHE[] = @get_scratch!("data")
    mkpath(CACHE[])
    return append!(MMDS, MixedModelsDatasets.datasets())
end

include("datasets.jl")
include("tagpad.jl")
include("movielens.jl")

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
