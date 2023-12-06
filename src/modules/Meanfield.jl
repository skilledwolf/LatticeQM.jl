
module Meanfield

    using SharedArrays, SparseArrays
    using LinearAlgebra

    using ..Utils
    import ..Spectrum: spectrum, chemicalpotential
    import ..TightBinding

    include("Meanfield/initialguess.jl")
    # export initialguess

    include("Meanfield/hartreefock.jl")
    # export hartreefock

    include("Meanfield/fixedpoint.jl")
    include("Meanfield/selfconsistent.jl")
    # export solveselfconsistent

end