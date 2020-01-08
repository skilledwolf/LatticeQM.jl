# __precompile__()
module LatticeQM # renamed from LatticeQM

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    # Export our base modules
    export Structure, BlochTools, TightBinding, KPM #, Plotting

    # Export high-level modules
    export Geometries2D, Materials

    # Dummy modules (provide stuff that is shared between modules)
    include("modules/DummySave.jl")
    using .DummySave
    export save

    ### Base modules
    include("modules/Structure.jl")
    using .Structure

    include("modules/BlochTools.jl")
    using .BlochTools
    export get_bands
    export solve_selfconsistent, initial_guess

    include("modules/TightBinding.jl")
    using .TightBinding
    export DenseHops, SparseHops, Hops, AbstractHops, hopdim, addhops!, addhops
    export get_operator, get_projector, get_hops, get_bloch, set_filling!

    # include("modules/Plotting.jl") # this module is deprecated in favor of using "recipes" with the package Plots.jl
    # using .Plotting

    include("modules/KPM.jl")
    using .KPM

    ### Higher-level modules
    include("modules/Geometries.jl")
    import .Geometries2D

    include("modules/Materials.jl")
    using .Materials

end

# using .Structure, .BlochTools, .TightBinding, .KPM
# using .Geometries2D, .Materials
