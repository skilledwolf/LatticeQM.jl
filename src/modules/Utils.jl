"""
General functions and definitions that don't necessarily make sense in the Algebra module.
"""
# Note: This module must NOT depend on any other module of the package.
module Utils

    using LinearAlgebra

    include("Utils/DummySave.jl")
    export DummySave


    include("Utils/legacymacros.jl")
    export @legacyalias, @legacymoved, @legacyremoved

    include("Utils/scalar2vector.jl")
    include("Utils/fermidirac.jl")

    include("Utils/rotate.jl")
    export RotationMatrix

    include("Utils/grid.jl")
end