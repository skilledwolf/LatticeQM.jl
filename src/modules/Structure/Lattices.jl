
"""
    Lattices

Provides the type `Lattice` and methods that act on this struct, see `?Structure.Lattices.Lattice`.



"""
module Lattices

    include("Lattices/type.jl")
    export Lattice

    include("Lattices/path.jl")
    export kpath

    include("Lattices/properties.jl")
    export latticedim, spacedim, countorbitals#, extraspacedim
    # export hasdimension, assertdimension
    export basis, getA, getB, coordinates, positions, allpositions, extracoordinates#, setextracoordinates!
    # export filterindices, filterpositions, filtercoordinates, fractionalize, fractionalize!, foldfractional

    include("Lattices/methods.jl")
    # export foldcoordinates!, rotatebasis!, rotatecoordinates!, translate!, mirrorZ, mirrorZ!
    # export newdimension!, mergelattices, mergelattices!

    include("Lattices/supercell.jl")
    export superlattice#, crop2unitcell, crop2unitcell!, bistack

    include("Lattices/reducedim.jl")
    # export reduceto0D, reduceto1D

    include("Lattices/neighbors.jl")
    # export getneighbors, getneighborcells, commonneighbor

    include("Lattices/twist.jl")
    # export twist

    include("Lattices/foldcell.jl")
    # export foldcell!, foldBZ!, foldPC!

    include("Lattices/fillregion.jl")
    
end