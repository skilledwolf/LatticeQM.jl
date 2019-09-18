__precompile__()
module LatticeQM # renamed from LatticeQM

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    # Export our base modules
    export Structure, BlochTools, TightBinding, Plotting

    # Export high-level modules
    export Geometries2D, Materials

    # Dummy modules (provide stuff that is shared between modules)
    include("modules/DummySave.jl")
    using .DummySave
    export save

    ### Base modules
    include("modules/Structure.jl")
    using .Structure

    include("modules/TightBinding.jl")
    using .TightBinding

    include("modules/BlochTools.jl")
    using .BlochTools

    include("modules/Plotting.jl")
    using .Plotting


    ### Higher-level modules
    include("modules/Geometries.jl")
    import .Geometries2D

    include("modules/Materials.jl")
    using .Materials

end
