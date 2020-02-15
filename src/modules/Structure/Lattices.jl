module Lattices

    using Plots
    using LinearAlgebra
    using RecursiveArrayTools

    using ...Utils
    using ..Paths

    include("Lattices/type.jl")
    export Lattice

    include("Lattices/properties.jl")
    export latticedim, extraspacedim, spacedim, countorbitals
    export hasdimension, assertdimension
    export getA, getB, coordinates, positions, allpositions, extrapositions, setextrapositions!
    export filterindices, filterpositions, filtercoordinates, fractionalize, fractionalize!, foldfractional

    include("Lattices/methods.jl")
    export foldcoordinates!, rotatebasis!, rotatecoordinates!, translate!, mirrorZ, mirrorZ!
    export newdimension!, mergelattices, mergelattices!

    include("Lattices/path.jl")
    export kpath


    include("Lattices/supercell.jl")
    export superlattice, turn0D!, repeat!, repeat, crop2unitcell, crop2unitcell!, bistack

    include("Lattices/plot.jl")
end