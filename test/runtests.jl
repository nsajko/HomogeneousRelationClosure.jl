using HomogeneousRelationClosure
using Test

function equal_vectors(l::AbstractVector, r::AbstractVector)
    if l != r
        throw("not equal")
    end
    l
end
function square_matrix_axis(a::AbstractMatrix)
    x = axes(a, 1)
    y = axes(a, 2)
    equal_vectors(x, y)
end
function square_matrix_axis(a::AbstractMatrix, bs::AbstractMatrix...)
    mapreduce(square_matrix_axis, equal_vectors, (a, bs...))
end
function relation_inclusion(l::AbstractMatrix, r::AbstractMatrix)
    square_matrix_axis(l, r)
    all(splat(≤), zip(l, r))
end
function relation_is_reflexive(a::AbstractMatrix)
    function f(i)
        Bool(a[i, i])::Bool
    end
    all(f, square_matrix_axis(a))
end
function relation_is_irreflexive(a::AbstractMatrix)
    function f(i)
        Bool(a[i, i])::Bool
    end
    all(!f, square_matrix_axis(a))
end
const relations = let
    function f_0(c)
        mat = Matrix{Float32}(undef, 0, 0)
        reshape(mat, :) .= c
        mat
    end
    function f_1(c)
        mat = Matrix{Float32}(undef, 1, 1)
        reshape(mat, :) .= c
        mat
    end
    function f_4(c)
        mat = Matrix{Float32}(undef, 2, 2)
        reshape(mat, :) .= c
        mat
    end
    function f_9(c)
        mat = Matrix{Float32}(undef, 3, 3)
        reshape(mat, :) .= c
        mat
    end
    i = Returns(0:1)
    relations_0 = Iterators.map(f_0, Iterators.product(ntuple(i, 0)...))
    relations_1 = Iterators.map(f_1, Iterators.product(ntuple(i, 1)...))
    relations_4 = Iterators.map(f_4, Iterators.product(ntuple(i, 4)...))
    relations_9 = Iterators.map(f_9, Iterators.product(ntuple(i, 9)...))
    Iterators.flatten((relations_0, relations_1, relations_4, relations_9))
end

@testset "HomogeneousRelationClosure.jl" begin
    @testset "reflexive" begin
        @testset "closure" begin
            for relation ∈ relations
                @test let a = copy(relation)
                    a === homogeneous_relation_reflexive_closure!(a)
                end
                @test let a = copy(relation)
                    relation_is_reflexive(homogeneous_relation_reflexive_closure!(a))
                end
                @test let a = copy(relation), b = copy(relation)
                    relation_inclusion(a, homogeneous_relation_reflexive_closure!(b))
                end
                @test let a = copy(relation)
                    homogeneous_relation_reflexive_closure!(a)
                    function g(x::AbstractMatrix)
                        (axes(relation) == axes(x)) && relation_is_reflexive(x) && relation_inclusion(relation, x)
                    end
                    all(Base.Fix1(relation_inclusion, a), Iterators.filter(g, relations))
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_reflexive_closure!([0 1])
        end
        @testset "reduction" begin
            for relation ∈ relations
                @test let a = copy(relation)
                    a === homogeneous_relation_reflexive_reduction!(a)
                end
                @test let a = copy(relation)
                    relation_is_irreflexive(homogeneous_relation_reflexive_reduction!(a))
                end
                @test let a = copy(relation), b = copy(relation)
                    relation_inclusion(homogeneous_relation_reflexive_reduction!(b), a)
                end
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_reflexive_closure!(a) ==
                    homogeneous_relation_reflexive_closure!(homogeneous_relation_reflexive_reduction!(b))
                end
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_reflexive_closure!(a)
                    homogeneous_relation_reflexive_reduction!(b)
                    function g(x::AbstractMatrix)
                        (axes(relation) == axes(x)) && (a == homogeneous_relation_reflexive_closure!(copy(x)))
                    end
                    all(Base.Fix1(relation_inclusion, b), Iterators.filter(g, relations))
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_reflexive_reduction!([0 1])
        end
    end
end

using Aqua: Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HomogeneousRelationClosure)
end
