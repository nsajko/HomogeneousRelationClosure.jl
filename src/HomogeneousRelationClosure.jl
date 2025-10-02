module HomogeneousRelationClosure
    export
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
    function homogeneous_relation_reflexive_closure!(a::AbstractMatrix)
        mutate_diagonal!(one, a)
    end
    function homogeneous_relation_reflexive_reduction!(a::AbstractMatrix)
        mutate_diagonal!(zero, a)
    end
    # TODO: transitive
    # TODO: symmetric
end
