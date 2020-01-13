# __precompile__()
module LatticeQM # renamed from LatticeQM

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    # Dummy modules (provide stuff that is shared between modules)
    include("modules/DummySave.jl")
    using .DummySave
    export save

    include("modules/Utils.jl")
    using .Utils

    include("modules/Algebra.jl")
    using .Algebra

    ### Base modules
    include("modules/Structure.jl")
    using .Structure
    export Structure

    include("modules/Spectrum.jl")
    using .Spectrum
    export Spectrum
    export get_bands

    include("modules/TightBinding.jl")
    using .TightBinding
    export TightBinding
    export DenseHops, SparseHops, Hops, AbstractHops, hopdim, addhops!, addhops
    export getoperator, getprojector, get_hops, get_bloch, set_filling!
    export initial_guess

    include("modules/Green.jl")
    using .Green
    export Green
    export density, density!
    export densitymatrix, densitymatrix!

    include("modules/Meanfield.jl")
    using .Meanfield
    export self_consistent # backwards compatibility
    export selfconsistent, hartreefock

    include("modules/KPM.jl")
    using .KPM
    export KPM

    ### Higher-level modules
    include("modules/Geometries.jl")
    import .Geometries2D
    export Geometries2D

    include("modules/Materials.jl")
    using .Materials
    export Materials

end

# using .Structure, .BlochTools, .TightBinding, .KPM
# using .Geometries2D, .Materials
