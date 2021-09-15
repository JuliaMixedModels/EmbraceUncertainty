module EU

    using Reexport
    @reexport using Books
    @reexport using CairoMakie
    @reexport using Chain
    @reexport using DataFrameMacros
    @reexport using DataFrames
    @reexport using MixedModels
    @reexport using MixedModelsMakie
    @reexport using PooledArrays
    @reexport using Statistics

    include("utilities.jl")
    include("intro.jl")
end
