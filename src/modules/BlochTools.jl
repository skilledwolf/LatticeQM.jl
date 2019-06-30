
module BlochTools

    using Base
    using ElasticArrays, SparseArrays
    using LinearAlgebra
    using Arpack

    using ..KSpace

    include("../misc/paulimatrices.jl")

    include("BlochTools/eigen.jl")
    include("BlochTools/density.jl")
    include("BlochTools/green.jl")

    include("BlochTools/meanfield.jl")

end
