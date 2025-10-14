using HomogeneousRelationClosure
using Test
using FixedSizeArrays: FixedSizeMatrix

module ExampleAlgebraicStructureTypes
    struct ExampleAlgebraicStructureType
        x::Float32
        global function create(x)
            new(x)
        end
    end
    const example_algebraic_structure = let
        function z(::ExampleAlgebraicStructureType)
            create(0)
        end
        function o(::ExampleAlgebraicStructureType)
            create(1)
        end
        function p(l::ExampleAlgebraicStructureType, r::ExampleAlgebraicStructureType)
            create(max(l.x, r.x))
        end
        function t(l::ExampleAlgebraicStructureType, r::ExampleAlgebraicStructureType)
            create(min(l.x, r.x))
        end
        function n(x::ExampleAlgebraicStructureType)
            create(1 - x.x)
        end
        (; zero = z, one = o, + = p, * = t, ! = n)
    end
    function Base.Bool(x::ExampleAlgebraicStructureType)  # used in tests
        Bool(x.x)
    end
    function Base.:(<)(l::ExampleAlgebraicStructureType, r::ExampleAlgebraicStructureType)  # used in tests
        l.x < r.x
    end
    function Base.transpose(x::ExampleAlgebraicStructureType)
        x
    end
end

module Booles
    export Boole
    """
        Boole

    Simplest non-trivial Boolean algebra, the two-element Boolean algebra.

    Subtypes `Integer`.

    Different from `Bool` in that `one(Boole) + one(Boole)` equals one instead of two.
    """
    struct Boole <: Integer
        b::Bool
        function Boole(b::Bool)
            new(b)
        end
    end
    function Boole(b::Int)
        Boole(Bool(b))
    end
    function Base.Bool(b::Boole)
        b.b
    end
    function Base.Int(b::Boole)
        Int(Bool(b))
    end
    function Base.:(+)(l::Boole, r::Boole)
        Boole(l.b | r.b)
    end
    function Base.:(*)(l::Boole, r::Boole)
        Boole(l.b & r.b)
    end
    function Base.:(<)(l::Boole, r::Boole)
        l.b < r.b
    end
    function Base.typemin(::Type{Boole})
        zero(Boole)
    end
    function Base.typemax(::Type{Boole})
        one(Boole)
    end
    function Base.:(!)(b::Boole)
        Boole(!(b.b))
    end
    function Base.show(io::IO, b::Boole)
        show(io, Boole)
        print(io, '(')
        show(io, Int(b))
        print(io, ')')
    end
end

using .Booles

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
const LogicalMatrix = Matrix{Boole}
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
    operations = let
        ops = [
            homogeneous_relation_reflexive_closure,
            homogeneous_relation_reflexive_closure!,
            homogeneous_relation_reflexive_reduction,
            homogeneous_relation_reflexive_reduction!,
            homogeneous_relation_symmetric_closure,
            homogeneous_relation_symmetric_closure!,
            homogeneous_relation_symmetric_reduction,
            homogeneous_relation_symmetric_reduction!,
            homogeneous_relation_transitive_closure,
            homogeneous_relation_transitive_closure!,
            homogeneous_relation_transitive_reduction_of_acyclic,
            homogeneous_relation_transitive_reduction_of_acyclic!,
        ]
        reshape(ops, (2, 2, :))
    end
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
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_reflexive_closure(a) == homogeneous_relation_reflexive_closure!(b)
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
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_reflexive_reduction(a) == homogeneous_relation_reflexive_reduction!(b)
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
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_symmetric_closure(a) == homogeneous_relation_symmetric_closure!(b)
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
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_symmetric_reduction(a) == homogeneous_relation_symmetric_reduction!(b)
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
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_transitive_closure(a) == homogeneous_relation_transitive_closure!(b)
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
                @test let a = copy(relation), b = copy(relation)
                    homogeneous_relation_transitive_reduction_of_acyclic(a) == homogeneous_relation_transitive_reduction_of_acyclic!(b)
                end
            end
            @test_throws DimensionMismatch homogeneous_relation_transitive_reduction_of_acyclic!([0 1])
        end
    end
    @testset "relation with custom scalars" begin
        relation = let ret = Matrix{ExampleAlgebraicStructureTypes.ExampleAlgebraicStructureType}(undef, 4, 4)
            z = ExampleAlgebraicStructureTypes.create(0)
            o = ExampleAlgebraicStructureTypes.create(1)
            ret .= Ref(z)
            ret[1, 2] = o
            ret[1, 3] = o
            ret[2, 4] = o
            ret[3, 4] = o
            ret
        end
        s = ExampleAlgebraicStructureTypes.example_algebraic_structure
        for operation ∈ operations
            @test ((@inferred operation(copy(relation), s)); true;)
        end
    end
    @testset "`FixedSizeMatrix`" begin
        relation = FixedSizeMatrix{Boole}(undef, 3, 3)
        for operation ∈ operations
            @test ((@inferred operation(copy(relation))); true;)
            @test operation(copy(relation)) isa FixedSizeMatrix{Boole}
        end
    end
end

using Aqua: Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HomogeneousRelationClosure)
end
