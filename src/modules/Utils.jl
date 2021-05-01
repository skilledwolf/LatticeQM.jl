"""
General functions and definitions that don't necessarily make sense anywhere else.
"""
# Note: This module must NOT depend on any other module of the package.
module Utils

    using LinearAlgebra

    # include("Utils/DummySave.jl")
    # export DummySave

    # include("Utils/legacymacros.jl")
    # export @legacyalias, @legacymoved, @legacyremoved

    include("Utils/pycall.jl")

    include("Utils/scalar2vector.jl")

    include("Utils/fermidirac.jl")
    include("Utils/paulimatrices.jl")
    export σ0, σ1, σ2, σ3, σX, σY, σZ, σUP, σDOWN, σPLUS, σMINUS, σs, spinorrotation


    """
        padvec(v::AbstractVector, d::Int)

    Make sure Vector v has length d, pad with zeros if needed.
    """
    function padvec(v::AbstractVector, d::Int)
        L = length(v)
        if L > d
            error("Vector exceeds specified length $d.")
        elseif L==d
            return v
        end
        return vcat(v,zeros(d-L))
    end
end
