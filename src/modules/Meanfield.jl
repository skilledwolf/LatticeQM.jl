
module Meanfield

    using ..Utils
    using ..Spectrum: spectrum, chemicalpotential
    using ..TightBinding
    using ..Green: densitymatrix!

    using SharedArrays

    include("Meanfield/types.jl")

    include("Meanfield/interactions.jl")
    export gethubbard, getcappedyukawa

    include("Meanfield/initialguess.jl")
    export initialguess

    include("Meanfield/fixedpoint.jl")
    include("Meanfield/selfconsistent.jl")
    export solveselfconsistent

    include("Meanfield/hartreefock.jl")
    export hartreefock, hartreefock_k



end