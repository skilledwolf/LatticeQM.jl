__precompile__()

module BlochTools

    using Base
    using ElasticArrays, SparseArrays
    using LinearAlgebra
    using Arpack

    using ..KSpace

    using ..TightBinding: get_bloch, get_operator

    export get_bands, spectrum, eigen, chemical_potential
    export solve_selfconsistent, initialize_ρ

    include("../misc/paulimatrices.jl")

    include("BlochTools/misc.jl")
    include("BlochTools/eigen.jl")
    include("BlochTools/density.jl")
    include("BlochTools/dos.jl")
    include("BlochTools/ldos.jl")
    include("BlochTools/density_matrix.jl")

    include("BlochTools/meanfield.jl")

end
