# __precompile__()

module Spectrum

    using Base
    using ElasticArrays, SparseArrays
    using LinearAlgebra
    using Arpack
    using Distributed

    # using ..KSpace

    using ..Utils
    using ..Algebra: σ0, σ1, σ2, σ3, σs
    using ..Structure.Paths: kIterable, DiscretePath, eachpoint, points, sumk

    export get_bands, spectrum, eigen, chemicalpotential
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
