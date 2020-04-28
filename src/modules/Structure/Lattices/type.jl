using ..Paths

"""
    Lattice

Type that contains all information about a lattice, such as lattice vectors
in ``Lattice.a``, fractional coordinates in ``Lattice.coordinates`` and
extra dimensions in ``Lattice.extradimensions``.

# Fields
- `a::Matrix{Float64}` :
  :math:`d \times d` matrix of lattice vectors a[:,i]
- `coordinates::Matrix{Float64}` :
  :math:`D \times d` matrix of fractional coordinates x[1:d,i]
  The coordinates x[d+1:D, i] are additional coordinates
- `extradimensions::Dict{String,Int}` :
  Keys are the names of the additional coordinates, e.g. "sublattice" or "z".
  The values are indices for the corresponding row in ``Lattice.coordinates``.
"""
mutable struct Lattice
    latticevectors::Matrix{Float64} # D × D
    orbitalcoordinates::Matrix{Float64} # d × N # fractional coordinates w.r.t. A
    extrapositions::Matrix{Float64} # D × N
    extradimensions::Dict{String,Int}
    specialpoints::LabeledPoints
end

kdict = LabeledPoints(Dict{String,AbstractVector{Float64}}(), Dict{String,String}(), Vector{String}([]))

Lattice(latticevectors::Matrix{Float64}) = Lattice(latticevectors, zeros(Float64, size(latticevectors,1), 1))

function Lattice(latticevectors::Matrix{Float64}, orbitalcoordinates::Matrix{Float64}; extradimensions::Vector{String}=Vector{String}(), specialpoints=kdict)
    @assert size(orbitalcoordinates,1) == size(latticevectors,2)
    # @assert size(latticevectors,1) == size(latticevectors,2)

    extrapositions = zeros(Float64, size(extradimensions), size(orbitalcoordinates,2)) # by default there is no perpendicular space
    Lattice(latticevectors, orbitalcoordinates, extrapositions; extradimensions=extradimensions_dict, specialpoints=specialpoints)
end

function Lattice(latticevectors::Matrix{Float64}, orbitalcoordinates::Matrix{Float64}, extrapositions::Matrix{Float64}; extradimensions::Vector{String}=Vector{String}(), specialpoints=kdict)
    @assert size(orbitalcoordinates,1) == size(latticevectors,2)
    # @assert size(latticevectors,1) == size(latticevectors,2)
    @assert size(extrapositions,2) == size(orbitalcoordinates,2)

    extradimensions_dict = Dict{String,UInt}(key=>index for (index,key) in enumerate(extradimensions))
    Lattice(latticevectors, orbitalcoordinates, extrapositions, extradimensions_dict, specialpoints)
end
