using HomogeneousRelationClosure
using Test
using Aqua

@testset "HomogeneousRelationClosure.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(HomogeneousRelationClosure)
    end
    # Write your tests here.
end
