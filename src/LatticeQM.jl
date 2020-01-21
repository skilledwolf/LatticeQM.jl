# __precompile__(false)
module LatticeQM

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    include("modules/Utils.jl")
    using .Utils
    export save

    include("modules/Algebra.jl")
    using .Algebra

    ### Base modules
    include("modules/Structure.jl")
    using .Structure
    export Structure
    export kpath

    include("modules/TightBinding.jl")
    using .TightBinding
    export TightBinding
    export DenseHops, SparseHops, Hops, AbstractHops, hopdim, addhops!, addhops
    export getoperator, getprojector, set_filling!
    export gethops, getbloch

    export get_hops, get_bloch # backwards compatibility
    export initial_guess # backwards compatibility

    include("modules/Spectrum.jl")
    using .Spectrum
    export Spectrum
    export getbands
    export get_bands # backwards compatibility

    include("modules/Green.jl")
    using .Green
    export Green
    export density, density!
    export densitymatrix, densitymatrix!

    include("modules/Operators.jl")
    using .Operators
    export Operators
    export getoperator, getprojector
    export trace, expval

    include("modules/LinearResponse.jl")
    using .LinearResponse

    include("modules/Meanfield.jl")
    using .Meanfield
    export self_consistent # backwards compatibility
    export selfconsistent, hartreefock
    export initialguess
    export gethubbard

    include("modules/KPM.jl")
    using .KPM
    export KPM

    ### Higher-level modules
    include("modules/Geometries.jl")
    import .Geometries2D
    export Geometries2D

    include("precompile.jl")
    _precompile_()
end
