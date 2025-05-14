
"""
    tagpad(v::AbstractVector{<:Integer}, ndig::Integer, tag::String="S"; pool::Bool=true)
    tagpad(v::AbstractVector{<:Integer}, tag::String="S"; pool::Bool=true)

Convert `v` to a vector of strings prepended with `tag` and padded to a constant string length of `ndig`,
which is evaluated as `maximum(ndigits, v)`, if not provided.

If `pool` is `true`, the default, the vector of strings is returned as a `PooledArray`.

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
