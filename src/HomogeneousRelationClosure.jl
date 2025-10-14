module HomogeneousRelationClosure
    module IndexArrays  # TODO: new package?
        export IndexArray
        mutable struct IndexArray{
            Index <: Integer,
            N,
            Wrapped <: (AbstractArray{T, N} where {T}),
        } <: AbstractArray{NTuple{N, Index}, N}
            const wrapped::Wrapped
            function IndexArray{Index}(wrapped::AbstractArray) where {Index <: Integer}
                n = ndims(wrapped)::Int
                if n < 0
                    throw(ArgumentError("negative ndims"))
                end
                Wrapped = typeof(wrapped)
                new{Index, n, Wrapped}(wrapped)
            end
        end
        function wrapped(a::IndexArray)
            getfield(a, 1)
        end
        function Base.getindex(::IndexArray{Index}, i::Vararg{Integer, N}) where {N, Index <: Integer}
            # TODO: check bounds?
            map(Index, i)
        end
        function Base.isassigned((@nospecialize a::IndexArray), (@nospecialize i::Integer))
            # TODO: check bounds?
            true
        end
        function Base.size(a::IndexArray)
            size(wrapped(a))
        end
        function Base.axes(a::IndexArray)
            axes(wrapped(a))
        end
        function Base.similar(a::IndexArray, ::Type{Element}, size::Tuple{Vararg{Int}}) where {Element}
            similar(wrapped(a), Element, size)
        end
        function Base.BroadcastStyle(
            ::Type{IndexArray{Index, N, Wrapped}},
        ) where {Index <: Integer, N, Wrapped <: (AbstractArray{T, N} where {T})}
            Base.BroadcastStyle(Wrapped)
        end
        function Base.copy(a::IndexArray)
            a
        end
        function Base.getindex(a::IndexArray, ::Colon)
            a
        end
    end
    export
        homogeneous_relation_transitive_closure,
        homogeneous_relation_transitive_reduction_of_acyclic,
        homogeneous_relation_symmetric_closure,
        homogeneous_relation_symmetric_reduction,
        homogeneous_relation_reflexive_closure,
        homogeneous_relation_reflexive_reduction,
        homogeneous_relation_transitive_closure!,
        homogeneous_relation_transitive_reduction_of_acyclic!,
        homogeneous_relation_symmetric_closure!,
        homogeneous_relation_symmetric_reduction!,
        homogeneous_relation_reflexive_closure!,
        homogeneous_relation_reflexive_reduction!
    using .IndexArrays
    const default_algebraic_structure = (; zero, one, +, *, !)
    function square_matrix_axis(a::AbstractMatrix)
        x = axes(a, 1)
        y = axes(a, 2)
        if x != y
            throw(DimensionMismatch())
        end
        x
    end
    function matrix_by_indices(
        predicate::Predicate,
        a::AbstractMatrix,
        ij::(NTuple{2,T} where {T <: Integer}),
        algebraic_structure = default_algebraic_structure,
    ) where {Predicate}
        (; zero, one) = algebraic_structure
        x = a[ij...]
        if predicate(ij...)
            one(x)
        else
            zero(x)
        end
    end
    function dot_product(
        l::AbstractVector,
        r::AbstractVector,
        algebraic_structure = default_algebraic_structure,
    )
        (; +, *, zero) = algebraic_structure
        if isempty(eachindex(l, r))
            zero(eltype(l)) + zero(eltype(r))
        else
            mapreduce(*, +, l, r)
        end
    end
    function square_matrix_product(
        l::AbstractMatrix,
        r::AbstractMatrix,
        algebraic_structure = default_algebraic_structure,
    )
        function func(ij::(NTuple{2,T} where {T <: Integer}))
            (i, j) = ij
            a = l[i, :]
            b = r[:, j]
            dot_product(a, b, algebraic_structure)
        end
        eachindex(l, r)
        index_matrix = IndexArray{Int}(l)
        func.(index_matrix)
    end
    """
        homogeneous_relation_reflexive_closure(a::AbstractMatrix, [algebraic_structure])

    The smallest reflexive relation that includes the input relation `a`.
    """
    function homogeneous_relation_reflexive_closure(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        square_matrix_axis(a)
        (; +, zero, one) = algebraic_structure
        function func(ij)  # identity
            matrix_by_indices(==, a, ij, algebraic_structure)
        end
        index_matrix = IndexArray{Int}(a)
        a .+ func.(index_matrix)
    end
    """
        homogeneous_relation_reflexive_reduction(a::AbstractMatrix, [algebraic_structure])

    The smallest relation included by the input relation `a` having the same reflexive closure as `a`.
    """
    function homogeneous_relation_reflexive_reduction(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        square_matrix_axis(a)
        (; *, zero, one) = algebraic_structure
        function func(ij)  # complement of identity
            matrix_by_indices(!=, a, ij, algebraic_structure)
        end
        index_matrix = IndexArray{Int}(a)
        a .* func.(index_matrix)
    end
    """
        homogeneous_relation_symmetric_closure(a::AbstractMatrix, [algebraic_structure])

    The smallest symmetric relation that includes `a`.
    """
    function homogeneous_relation_symmetric_closure(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        square_matrix_axis(a)
        (; +) = algebraic_structure
        c = transpose(a)  # converse relation
        a .+ c
    end
    """
        homogeneous_relation_symmetric_reduction(a::AbstractMatrix, [algebraic_structure])

    Smallest relation included by the input relation `a` having the same symmetric closure as `a`.
    """
    function homogeneous_relation_symmetric_reduction(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        square_matrix_axis(a)
        (; +, *, !, zero, one) = algebraic_structure
        function func(ij)  # upper triangular
            matrix_by_indices(<=, a, ij, algebraic_structure)
        end
        index_matrix = IndexArray{Int}(a)
        c = transpose(a)  # converse relation
        a .* (func.(index_matrix) .+ (!).(a .* c))
    end
    """
        homogeneous_relation_transitive_closure(a::AbstractMatrix, [algebraic_structure])

    The smallest transitive relation that includes the input relation `a`.
    """
    function homogeneous_relation_transitive_closure(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        axis = square_matrix_axis(a)
        (; +) = algebraic_structure
        ret = a
        matrix_power = a
        for _ ∈ axis[2:end]
            matrix_power = square_matrix_product(matrix_power, a, algebraic_structure)
            ret = ret .+ matrix_power
        end
        ret
    end
    """
        homogeneous_relation_transitive_reduction_of_acyclic(a::AbstractMatrix, [algebraic_structure])

    The smallest relation included by the input relation `a` having the same transitive closure as `a`.

    Silently produces incorrect results in the presence of cycles.
    """
    function homogeneous_relation_transitive_reduction_of_acyclic(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        square_matrix_axis(a)
        (; *, !) = algebraic_structure
        # TODO: optimize: eliminate one matrix multiplication
        b = homogeneous_relation_transitive_closure(a, algebraic_structure)
        ab = square_matrix_product(a, b, algebraic_structure)
        a .* (!).(ab)
    end
    function mutate_diagonal!(f::F, a::AbstractMatrix) where {F}
        axis = square_matrix_axis(a)
        for i ∈ axis
            x = a[i, i]
            y = f(x)
            a[i, i] = y
        end
        a
    end
    """
        homogeneous_relation_reflexive_closure!(a::AbstractMatrix, [algebraic_structure])

    The smallest reflexive relation that includes the input relation `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_reflexive_closure!(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        (; one) = algebraic_structure
        mutate_diagonal!(one, a)
    end
    """
        homogeneous_relation_reflexive_reduction!(a::AbstractMatrix, [algebraic_structure])

    The smallest relation included by the input relation `a` having the same reflexive closure as `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_reflexive_reduction!(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        (; zero) = algebraic_structure
        mutate_diagonal!(zero, a)
    end
    """
        homogeneous_relation_symmetric_closure!(a::AbstractMatrix, [algebraic_structure])

    The smallest symmetric relation that includes the input relation `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_symmetric_closure!(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        axis = square_matrix_axis(a)
        for I ∈ eachindex(axis)
            i = axis[I]
            for J ∈ (I + 1):lastindex(axis)
                j = axis[J]
                x_ij = a[i, j]
                x_ji = a[j, i]
                y = let (; +) = algebraic_structure
                    x_ij + x_ji
                end
                a[i, j] = y
                a[j, i] = y
            end
        end
        a
    end
    """
        homogeneous_relation_symmetric_reduction!(a::AbstractMatrix, [algebraic_structure])

    Smallest relation included by the input relation `a` having the same symmetric closure as `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_symmetric_reduction!(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        axis = square_matrix_axis(a)
        for I ∈ eachindex(axis)
            i = axis[I]
            for J ∈ (I + 1):lastindex(axis)
                j = axis[J]
                x_ij = a[i, j]
                x_ji = a[j, i]
                y_ji = let (; *, !) = algebraic_structure
                    (!x_ij) * x_ji
                end
                a[j, i] = y_ji
            end
        end
        a
    end
    """
        homogeneous_relation_transitive_closure!(a::AbstractMatrix, [algebraic_structure])

    The smallest transitive relation that includes the input relation `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_transitive_closure!(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        axis = square_matrix_axis(a)
        for k ∈ axis
            for i ∈ axis
                for j ∈ axis
                    x_ij = a[i, j]
                    x_ik = a[i, k]
                    x_kj = a[k, j]
                    y_ij = let (; +, *) = algebraic_structure
                        x_ij + (x_ik * x_kj)
                    end
                    a[i, j] = y_ij
                end
            end
        end
        a
    end
    """
        homogeneous_relation_transitive_reduction_of_acyclic!(a::AbstractMatrix, [algebraic_structure])

    The smallest relation included by the input relation `a` having the same transitive closure as `a`.

    Silently produces incorrect results in the presence of cycles.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_transitive_reduction_of_acyclic!(a::AbstractMatrix, algebraic_structure = default_algebraic_structure)
        homogeneous_relation_transitive_closure!(a, algebraic_structure)
        axis = square_matrix_axis(a)
        rax = reverse(axis)
        for k ∈ rax
            for i ∈ axis
                for j ∈ axis
                    x_ij = a[i, j]
                    x_ik = a[i, k]
                    x_kj = a[k, j]
                    y_ij = let (; *, !) = algebraic_structure
                        x_ij * !(x_ik * x_kj)
                    end
                    a[i, j] = y_ij
                end
            end
        end
        a
    end
end
