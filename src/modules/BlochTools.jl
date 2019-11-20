# __precompile__()

module BlochTools

    using Base
    using ElasticArrays, SparseArrays
    using LinearAlgebra
    using Arpack

    # using ..KSpace

    using ..Structure: kIterable, DiscretePath, eachpoint, points, sumk

    export get_bands, spectrum, eigen, chemical_potential
    export solve_selfconsistent

    include("../misc/paulimatrices.jl")

    include("BlochTools/types.jl")

    include("BlochTools/misc.jl")
    include("BlochTools/eigen.jl")
    include("BlochTools/density.jl")
    include("BlochTools/dos.jl")
    include("BlochTools/ldos.jl")
    include("BlochTools/density_matrix.jl")
    include("BlochTools/berry.jl")

    include("BlochTools/meanfield.jl")

end
