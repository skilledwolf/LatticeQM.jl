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

    include("rotate.jl")

    include("Paths.jl")
    import .Paths

    include("grid.jl")

    include("Lattices.jl")
    import .Lattices
    import .Lattices: Lattice, kpath
    export Lattice, kpath

    include("Geometries.jl")
    import .Geometries

################################################################################
################################################################################
end
