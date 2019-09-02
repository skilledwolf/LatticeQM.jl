__precompile__()
module LatticeQM # renamed from LatticeQM

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    # Export our base modules
    export Structure, KSpace, BlochTools, TightBinding, Plotting

    ### Base modules
    include("modules/Structure.jl")
    using .Structure

    include("modules/KSpace.jl")
    using .KSpace

    include("modules/TightBinding.jl")
    using .TightBinding

    include("modules/BlochTools.jl")
    using .BlochTools

    include("modules/Plotting.jl")
    using .Plotting


    ### Higher-level modules
    include("modules/Geometries.jl")

end
