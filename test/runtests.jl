using HomogeneousRelationClosure
using Test

@testset "HomogeneousRelationClosure.jl" begin
    # Write your tests here.
end

using Aqua: Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HomogeneousRelationClosure)
end
