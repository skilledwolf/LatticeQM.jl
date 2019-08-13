__precompile__()

module BlochTools

    using Base
    using ElasticArrays, SparseArrays
    using LinearAlgebra
    using Arpack

    using ..KSpace

    using ..TightBinding: build_BlochH, get_operator

    include("../misc/paulimatrices.jl")

    include("BlochTools/misc.jl")
    include("BlochTools/eigen.jl")
    include("BlochTools/density.jl")
    include("BlochTools/dos.jl")
    include("BlochTools/ldos.jl")
    include("BlochTools/density_matrix.jl")

    include("BlochTools/meanfield.jl")

end
