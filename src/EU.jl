module EU

import Books

using Reexport
@reexport using Books
@reexport using CairoMakie
@reexport using CategoricalArrays
@reexport using Chain
@reexport using DataFrameMacros
@reexport using DataFrames
@reexport using MixedModels
@reexport using MixedModelsMakie
@reexport using PooledArrays
@reexport using Statistics

include("utilities.jl")
include("intro.jl")

"""
    build()

This method is called during CI.
"""
function build()
    println("Building EU")
    fail_on_error = true
    Books.gen(; fail_on_error)
    Books.build_all(; fail_on_error)
end

end # module
