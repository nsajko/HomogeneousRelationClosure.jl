using HomogeneousRelationClosure
using Test
using BooleanSemiring: B

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
function reflexivity_predicate(a::AbstractMatrix, b::Bool)
    function f(i)
        Bool(a[i, i])::Bool
    end
    all(Base.Fix1(==, b) ∘ f, square_matrix_axis(a))
end
function relation_is_reflexive(a::AbstractMatrix)
    reflexivity_predicate(a, true)
end
function relation_is_irreflexive(a::AbstractMatrix)
    reflexivity_predicate(a, false)
end
function relation_is_symmetric(a::AbstractMatrix)
    a == transpose(a)
end
function relation_is_antisymmetric(a::AbstractMatrix)
    function f((i, j))
        Bool(a[i, j] * a[j, i])::Bool
    end
    axis = square_matrix_axis(a)
    all(!f, Iterators.filter(splat(!=), Iterators.product(axis, axis)))
end
function transitivity_predicate(a::AbstractMatrix, b::Bool)
    function f((i, _, k))
        Bool(a[i, k])::Bool
    end
    function g((i, j, k))
        Bool(a[i, j] * a[j, k])::Bool
    end
    axis = square_matrix_axis(a)
    all(Base.Fix1(==, b) ∘ f, Iterators.filter(g, Iterators.product(axis, axis, axis)))
end
function relation_is_transitive(a::AbstractMatrix)
    transitivity_predicate(a, true)
end
function relation_is_antitransitive(a::AbstractMatrix)
    transitivity_predicate(a, false)
end
function relation_is_dag(a::AbstractMatrix)
    # Kahn's algorithm for topological sorting
    function f(i)
        all(iszero, @view a[:, i])
    end
    axis = square_matrix_axis(a)
    S = collect(eltype(axis), Iterators.filter(f, axis))
    b = copy(a)
    while !isempty(S)
        n = pop!(S)
        function g(m)
            isone(b[n, m])::Bool
        end
        for m ∈ Iterators.filter(g, axis)
            b[n, m] = zero(b[n, m])
            if all(iszero, @view b[:, m])
                push!(S, m)
            end
        end
    end
    iszero(b)
end
const LogicalMatrix = Matrix{B}
const relations = let
    function f(c)
        len = length(c)
        n = isqrt(len)
        if n * n != len
            throw(DimensionMismatch())
        end
        mat = LogicalMatrix(undef, n, n)
        reshape(mat, :) .= c
        mat
    end
    function g(n)
        Iterators.map(f, Iterators.product(ntuple(Returns(0:1), n)...))
    end
    Iterators.flatten(Iterators.map(g, (0, 1, 4, 9)))
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
                @test let a = copy(relation)
                    relation_inclusion(relation, homogeneous_relation_reflexive_closure!(a))
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
                @test let a = copy(relation)
                    relation_inclusion(homogeneous_relation_reflexive_reduction!(a), relation)
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
    @testset "symmetric" begin
        @testset "closure" begin
            for relation ∈ relations
                @test let a = copy(relation)
                    a === homogeneous_relation_symmetric_closure!(a)
                end
                @test let a = copy(relation)
                    relation_is_symmetric(homogeneous_relation_symmetric_closure!(a))
                end
                @test let a = copy(relation)
                    relation_inclusion(relation, homogeneous_relation_symmetric_closure!(a))
                end
                @test let a = copy(relation)
                    homogeneous_relation_symmetric_closure!(a)
                    function g(x::AbstractMatrix)
                        (axes(relation) == axes(x)) && relation_is_symmetric(x) && relation_inclusion(relation, x)
                    end
                    all(Base.Fix1(relation_inclusion, a), Iterators.filter(g, relations))
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_symmetric_closure!([0 1])
        end
        @testset "reduction" begin
            for relation ∈ relations
                @test let a = copy(relation)
                    a === homogeneous_relation_symmetric_reduction!(a)
                end
                @test let a = copy(relation)
                    relation_is_antisymmetric(homogeneous_relation_symmetric_reduction!(a))
                end
                @test let a = copy(relation)
                    relation_inclusion(homogeneous_relation_symmetric_reduction!(a), relation)
                end
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_symmetric_closure!(a) ==
                    homogeneous_relation_symmetric_closure!(homogeneous_relation_symmetric_reduction!(b))
                end
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_symmetric_closure!(a)
                    homogeneous_relation_symmetric_reduction!(b)
                    function f(x::AbstractMatrix)
                        count(isone, b) ≤ count(isone, x)
                    end
                    function g(x::AbstractMatrix)
                        (axes(relation) == axes(x)) && (a == homogeneous_relation_symmetric_closure!(copy(x)))
                    end
                    all(f, Iterators.filter(g, relations))
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_symmetric_reduction!([0 1])
        end
    end
    @testset "transitive" begin
        @testset "closure" begin
            for relation ∈ relations
                @test let a = copy(relation)
                    a === homogeneous_relation_transitive_closure!(a)
                end
                @test let a = copy(relation)
                    relation_is_transitive(homogeneous_relation_transitive_closure!(a))
                end
                @test let a = copy(relation)
                    relation_inclusion(relation, homogeneous_relation_transitive_closure!(a))
                end
                @test let a = copy(relation)
                    homogeneous_relation_transitive_closure!(a)
                    function g(x::AbstractMatrix)
                        (axes(relation) == axes(x)) && relation_is_transitive(x) && relation_inclusion(relation, x)
                    end
                    all(Base.Fix1(relation_inclusion, a), Iterators.filter(g, relations))
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_transitive_closure!([0 1])
        end
        @testset "reduction for acyclic" begin
            for relation ∈ Iterators.filter(relation_is_dag, relations)
                @test let a = copy(relation)
                    a === homogeneous_relation_transitive_reduction_of_acyclic!(a)
                end
                @test let a = copy(relation)
                    relation_is_antitransitive(homogeneous_relation_transitive_reduction_of_acyclic!(a))
                end
                @test let a = copy(relation)
                    relation_inclusion(homogeneous_relation_transitive_reduction_of_acyclic!(a), relation)
                end
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_transitive_closure!(a) ==
                    homogeneous_relation_transitive_closure!(homogeneous_relation_transitive_reduction_of_acyclic!(b))
                end
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_transitive_closure!(a)
                    homogeneous_relation_transitive_reduction_of_acyclic!(b)
                    function f(x::AbstractMatrix)
                        count(isone, b) ≤ count(isone, x)
                    end
                    function g(x::AbstractMatrix)
                        (axes(relation) == axes(x)) && (a == homogeneous_relation_transitive_closure!(copy(x)))
                    end
                    all(f, Iterators.filter(g, relations))
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_transitive_reduction_of_acyclic!([0 1])
        end
    end
end

using Aqua: Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HomogeneousRelationClosure)
end
