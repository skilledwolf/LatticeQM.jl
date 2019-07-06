# __precompile__()
module LatticeToolbox

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    # Export the base modules
    export LatticeTools, KSpace, BlochTools, TightBinding, Plotting

    ### Base modules
    include("modules/LatticeTools.jl")
    using .LatticeTools

    include("modules/KSpace.jl")
    using .KSpace

    include("modules/BlochTools.jl")
    using .BlochTools

    include("modules/TightBinding.jl")
    using .TightBinding

    include("modules/Plotting.jl")
    using .Plotting

end
