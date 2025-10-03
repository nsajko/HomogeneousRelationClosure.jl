module HomogeneousRelationClosure
    export
        homogeneous_relation_transitive_closure!,
        homogeneous_relation_transitive_reduction_of_acyclic!,
        homogeneous_relation_symmetric_closure!,
        homogeneous_relation_symmetric_reduction!,
        homogeneous_relation_reflexive_closure!,
        homogeneous_relation_reflexive_reduction!
    function square_matrix_axis(a::AbstractMatrix)
        x = axes(a, 1)
        y = axes(a, 2)
        if x != y
            throw(DimensionMismatch())
        end
        x
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
        homogeneous_relation_reflexive_closure!(a::AbstractMatrix)

    The smallest reflexive relation that includes `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_reflexive_closure!(a::AbstractMatrix)
        mutate_diagonal!(one, a)
    end
    """
        homogeneous_relation_reflexive_reduction!(a::AbstractMatrix)

    The smallest relation included by `a` that has the same reflexive closure as `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_reflexive_reduction!(a::AbstractMatrix)
        mutate_diagonal!(zero, a)
    end
    """
        homogeneous_relation_symmetric_closure!(a::AbstractMatrix)

    The smallest symmetric relation that includes `a`.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_symmetric_closure!(a::AbstractMatrix)
        axis = square_matrix_axis(a)
        for I ∈ eachindex(axis)
            i = axis[I]
            for J ∈ (I + 1):lastindex(axis)
                j = axis[J]
                x_ij = a[i, j]
                x_ji = a[j, i]
                y = x_ij + x_ji
                a[i, j] = y
                a[j, i] = y
            end
        end
        a
    end
    """
        `homogeneous_relation_symmetric_reduction!(a::AbstractMatrix)`

    Smallest relation included by the input relation that has the same symmetric closure as the input relation.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_symmetric_reduction!(a::AbstractMatrix)
        axis = square_matrix_axis(a)
        for I ∈ eachindex(axis)
            i = axis[I]
            for J ∈ (I + 1):lastindex(axis)
                j = axis[J]
                x_ij = a[i, j]
                x_ji = a[j, i]
                y_ji = (!x_ij) * x_ji
                a[j, i] = y_ji
            end
        end
        a
    end
    """
        homogeneous_relation_transitive_closure!(a::AbstractMatrix)

    The smallest transitive relation that includes the input relation.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_transitive_closure!(a::AbstractMatrix)
        axis = square_matrix_axis(a)
        for k ∈ axis
            for i ∈ axis
                for j ∈ axis
                    x_ik = a[i, k]
                    x_kj = a[k, j]
                    y_ij = x_ik * x_kj
                    a[i, j] += y_ij
                end
            end
        end
        a
    end
    """
        homogeneous_relation_transitive_reduction_of_acyclic!(a::AbstractMatrix)

    The smallest relation included by the input relation that has the same transitive closure as the input relation.

    Silently produces incorrect results in the presence of cycles.

    Mutate `a`. Return `a`.
    """
    function homogeneous_relation_transitive_reduction_of_acyclic!(a::AbstractMatrix)
        homogeneous_relation_transitive_closure!(a)
        axis = square_matrix_axis(a)
        rax = reverse(axis)
        for k ∈ rax
            for i ∈ axis
                for j ∈ axis
                    x_ik = a[i, k]
                    x_kj = a[k, j]
                    y_ij = x_ik * x_kj
                    a[i, j] *= !y_ij
                end
            end
        end
        a
    end
end
