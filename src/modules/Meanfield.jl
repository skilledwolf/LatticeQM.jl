
module Meanfield

    using SharedArrays, SparseArrays
    using LinearAlgebra

    using ..Utils
    using ..Algebra
    using ..Spectrum: spectrum, chemicalpotential
    using ..TightBinding
    using ..Green: densitymatrix!

    include("Meanfield/types.jl")

    include("Meanfield/interactions.jl")
    export gethubbard, getcappedyukawa

    include("Meanfield/initialguess.jl")
    export initialguess

    include("Meanfield/fixedpoint.jl")
    include("Meanfield/selfconsistent.jl")
    export solveselfconsistent, solvehartreefockmf_k

    include("Meanfield/hartreefock.jl")
    export hartreefock, hartreefock_k

end