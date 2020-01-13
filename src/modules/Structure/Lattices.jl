module Lattices

    using Plots
    using LinearAlgebra

    using ...Utils
    using ..Paths

    include("Lattices/type.jl")
    export Lattice

    include("Lattices/properties.jl")
    export latticedim, extraspacedim, spacedim, countatoms
    export getA, getB, coordinates, positions, allpositions, extrapositions, setextrapositions
    export filterindices, filterpositions, filtercoordinates, fractionalize, fractionalize!, foldfractional

    include("Lattices/methods.jl")
    export foldcoordinates!, rotatebasis!, rotatecoordinates!, translate!, mirrorZ, mirrorZ!
    export newdimension!, mergelattices, mergelattices!

    include("Lattices/path.jl")
    export kpath

    include("Lattices/plot.jl")
end