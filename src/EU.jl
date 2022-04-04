module EU

using Books: Books

using Reexport
@reexport using AlgebraOfGraphics
@reexport using Arrow
@reexport using Books
@reexport using CairoMakie
@reexport using CategoricalArrays
@reexport using Chain
@reexport using DataFrameMacros
@reexport using DataFrames
@reexport using LinearAlgebra
@reexport using MixedModels
@reexport using MixedModelsMakie
@reexport using PooledArrays
@reexport using Random
@reexport using StandardizedPredictors
@reexport using Statistics
@reexport using StatsBase

include("utilities.jl")
include("intro.jl")
include("multiple.jl")
include("longitudinal.jl")

"""
    build()

This method is called during CI.
"""
function build()
    println("Building EU")
    fail_on_error = true
    mkpath(Books.BUILD_DIR)
    Books.gen(; fail_on_error)
    # Books.html(; fail_on_error) should be replaced by
    # Books.build_all(; fail_on_error)
    # once the figures such as "Multiple-fm05Limage" have been replaced.
    # Verifies cross references.
    return Books.html(; fail_on_error)
end

end # module
