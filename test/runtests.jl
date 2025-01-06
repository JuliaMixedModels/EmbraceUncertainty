using Aqua
using EmbraceUncertainty
using Test
using TestSetExtensions

@testset ExtendedTestSet "EmbraceUncertainty" begin
    @testset "Aqua" begin
        # we can't test stale deps because the dependencies aren't 
        # just for the package, but also for entire book
        Aqua.test_all(EmbraceUncertainty; ambiguities=false, stale_deps=false)
    end
end
