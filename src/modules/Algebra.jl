"""
    Algebra

This module intended to be a collection of matrix methods (and consistent wrappers).
It is meant to be fully independent of the other modules of this package.
"""
module Algebra

    using LinearAlgebra
    using Arpack
    # using KrylovKit


    include("Algebra/paulimatrices.jl")
    export σ0, σ1, σ2, σ3, σX, σY, σZ, σs
    export spinorrotation

    include("Algebra/eigen.jl")

end