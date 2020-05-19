__precompile__(true)
module LatticeQM

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    using KPM # this is one of my packages

    include("modules/Utils.jl")
    using .Utils
    using .DummySave
    export save, savedlm

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
    export getoperator, getprojector
    export gethops, getbloch

    include("modules/Spectrum.jl")
    using .Spectrum
    export Spectrum
    export getbands

    include("modules/Green.jl")
    using .Green
    export Green
    export density, density!
    export densitymatrix, densitymatrix!

    include("modules/Operators.jl")
    using .Operators
    export Operators
    export getoperator, getprojector
    export setfilling!
    export trace, expval

    include("modules/LinearResponse.jl")
    using .LinearResponse

    include("modules/Meanfield.jl")
    using .Meanfield
    export selfconsistent, hartreefock
    export initialguess
    export gethubbard

#     # include("modules/KPM.jl") # was moved to a separate package
#     # using .KPM
#     # export KPM

    ### Higher-level modules
    include("modules/Geometries.jl")
    import .Geometries2D
    export Geometries2D

    ### Plotting recipes
    using RecipesBase
    include("plotting/Lattice.jl")
    include("plotting/Bands.jl")

# #     include("precompile.jl")
# #     _precompile_()
end
