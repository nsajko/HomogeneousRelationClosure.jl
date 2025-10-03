# HomogeneousRelationClosure

[![Build Status](https://github.com/nsajko/HomogeneousRelationClosure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nsajko/HomogeneousRelationClosure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Package version](https://juliahub.com/docs/General/HomogeneousRelationClosure/stable/version.svg)](https://juliahub.com/ui/Packages/General/HomogeneousRelationClosure)
[![Package dependencies](https://juliahub.com/docs/General/HomogeneousRelationClosure/stable/deps.svg)](https://juliahub.com/ui/Packages/General/HomogeneousRelationClosure?t=2)
[![Coverage](https://codecov.io/gh/nsajko/HomogeneousRelationClosure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nsajko/HomogeneousRelationClosure.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/H/HomogeneousRelationClosure.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/H/HomogeneousRelationClosure.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A software package for the Julia programming language implementing *closure* and *reduction* operators for some properties of homogeneous binary relations, such as *reflexivity*, *symmetry* or *transitivity*.

Each homogeneous binary relation is assumed to be represented by a (square, logical) adjacency matrix (`AbstractMatrix`), whose elements support multiplication/conjunction (`*`), addition/disjunction (`+`) and negation (`!`). The element type is intended to be a Boolean algebra, such as provided by the package [TwoElementBooleanAlgebra.jl](https://github.com/nsajko/TwoElementBooleanAlgebra.jl). A fuzzy logic implementing the same operations should also work, in place of a Boolean logic (no idea if there's any application for that).

## Functionality

* reflexivity

    * *reflexive closure*

        * The smallest reflexive relation that includes the input relation.

        * Unique.

        * `homogeneous_relation_reflexive_closure!(::AbstractMatrix)`

    * *reflexive reduction*

        * The smallest relation included by the input relation that has the same reflexive closure as the input relation.

        * Unique.

        * `homogeneous_relation_reflexive_reduction!(::AbstractMatrix)`

* symmetry

    * *symmetric closure*

        * The smallest symmetric relation that includes the input relation.

        * Unique.

        * `homogeneous_relation_symmetric_closure!(::AbstractMatrix)`

    * *symmetric reduction*

        * Smallest relation included by the input relation that has the same symmetric closure as the input relation.

        * Not unique.

        * `homogeneous_relation_symmetric_reduction!(::AbstractMatrix)`

* transitivity

    * *transitive closure*

        * The smallest transitive relation that includes the input relation.

        * Unique.

        * `homogeneous_relation_transitive_closure!(::AbstractMatrix)`

    * *transitive reduction*

        * Smallest relation included by the input relation that has the same transitive closure as the input relation.

        * Not unique, except for directed acyclic graphs (DAGs).

        * `homogeneous_relation_transitive_reduction_of_acyclic!(::AbstractMatrix)`: warning: silently produces incorrect results in the presence of cycles
