# __precompile__()

"""
    Structure 

Provides the struct `Lattices.Lattice` to define and manipulate lattices, 
and `Paths.DiscretePath` to deal with discretized paths.

Check out the submodules:
- Lattices
- Paths
- Geometries

"""
module Structure
################################################################################
################################################################################

    include("Structure/rotate.jl")

    include("Structure/Paths.jl")
    import .Paths

    include("Structure/grid.jl")

    include("Structure/Lattices.jl")
    import .Lattices
    import .Lattices: Lattice, kpath
    export Lattice, kpath

    include("Structure/Geometries.jl")
    import .Geometries

################################################################################
################################################################################
end
