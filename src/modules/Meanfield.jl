
module Meanfield

    using SharedArrays, SparseArrays
    using LinearAlgebra

    using ..Utils
    import ..Spectrum: spectrum, chemicalpotential
    import ..TightBinding

    include("Meanfield/initialguess.jl")
    export initialguess

    include("Meanfield/fixedpoint.jl")
    include("Meanfield/selfconsistent.jl")
    export solveselfconsistent, solvehartreefockmf_k

    include("Meanfield/hartreefock.jl")
    export hartreefock, hartreefock_k

end