# __precompile__()
module LatticeQM # renamed from LatticeToolbox

    using Base, ElasticArrays, SparseArrays, LinearAlgebra

    # Export the base modules
    export LatticeTools, KSpace, BlochTools, TightBinding, Plotting

    ### Base modules
    include("modules/LatticeTools.jl")
    using .LatticeTools

    include("modules/KSpace.jl")
    using .KSpace

    include("modules/TightBinding.jl")
    using .TightBinding

    include("modules/BlochTools.jl")
    using .BlochTools

    include("modules/Plotting.jl")
    using .Plotting

end
