# __precompile__()

module Spectrum

    using ElasticArrays, SparseArrays
    using LinearAlgebra
    using Arpack
    using Distributed
    using SharedArrays

    # using ..KSpace

    using ..Utils
    using ..Algebra: σ0, σ1, σ2, σ3, σs
    using ..Structure: Lattice
    using ..Structure.Paths: kIterable, DiscretePath, eachpoint, points, sumk

    export getbands, spectrum, eigen, chemicalpotential
    export chemical_potential # backwards compatibility

    include("Spectrum/types.jl")

    include("Spectrum/spectrum.jl")
    include("Spectrum/bands.jl")
    include("Spectrum/chemicalpotential.jl")
    include("Spectrum/bandgap.jl")
    include("Spectrum/groundstatesum.jl")

    include("Spectrum/dos.jl")
    include("Spectrum/ldos.jl")

    include("Spectrum/berry.jl")

    #include("BlochTools/meanfield.jl")
end
