"""
This is intended to be a collection of matrix methods (and wrappers).
It is meant to fully independent of the other modules of this package.
"""
module Algebra

    using LinearAlgebra

    include("Algebra/paulimatrices.jl")
    export σ0, σ1, σ2, σ3, σs

    include("Algebra/eigen.jl")

end