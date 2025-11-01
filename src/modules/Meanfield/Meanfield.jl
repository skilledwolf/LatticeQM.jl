
module Meanfield

    using SharedArrays, SparseArrays
    using LinearAlgebra

    import LatticeQM.Utils
    import LatticeQM.TightBinding

    include("initialguess.jl")
    # export initialguess

    include("hartreefock.jl")
    # export hartreefock

    include("fixedpoint.jl")
    include("selfconsistent.jl")
    include("selfconsistent_purification.jl")
    # export solveselfconsistent

end