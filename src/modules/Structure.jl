# __precompile__()

"""
    Structure 

Provides the struct `Lattices.Lattice` to define and manipulate lattices, 
and `Paths.DiscretePath` to deal with discretized paths.

Check out the submodules:
- Lattices
- Paths

"""
module Structure
################################################################################
################################################################################

    include("Structure/Paths.jl")
    using .Paths

    include("Structure/Lattices.jl")
    using .Lattices
    export Lattice
    export kpath

################################################################################
################################################################################
end
