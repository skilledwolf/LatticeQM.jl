
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
mutable struct Lattice{T<:AbstractMatrix{Float64}}

    A::T # d × d
    atoms::T # d × N # fractional coordinates w.r.t. A

    atoms_aux::T # D × N
    extradimensions::Dict{String,Int}

    highsymmetrypoints::LabeledPoints
end

kdict = LabeledPoints(Dict{String,AbstractVector{Float64}}(), Dict{String,String}(), Vector{String}([]))

Lattice(A::T) where {T<:AbstractMatrix{Float64}} = Lattice{T}(A, zeros(Float64, size(A,1), 1))

function Lattice(A::T, atoms::T; extradimensions::Vector{String}=Vector{String}(), highsymmetrypoints=kdict) where {T<:AbstractMatrix{Float64}}
    @assert size(atoms,1) == size(A,1)
    @assert size(A,1) == size(A,2)

    L = size(extradimensions)
    atoms_aux = zeros(Float64, L, size(atoms,2)) # by default there is no perpendicular space

    # extradimensions_dict = Dict{String,Int}(key=>index for (index,key) in enumerate(extradimensions))
    # Lattice{T}(A, atoms, atoms_aux, extradimensions_dict)

    Lattice(A, atoms, atoms_aux; extradimensions=extradimensions_dict, highsymmetrypoints=highsymmetrypoints)
end

function Lattice(A::T, atoms::T, atoms_aux::T; extradimensions::Vector{String}=Vector{String}(), highsymmetrypoints=kpoints) where {T<:AbstractMatrix{Float64}}
    @assert size(atoms,1) == size(A,1)
    @assert size(A,1) == size(A,2)
    @assert size(atoms_aux,2) == size(atoms,2)

    extradimensions_dict = Dict{String,UInt}(key=>index for (index,key) in enumerate(extradimensions))

    Lattice{T}(A, atoms, atoms_aux, extradimensions_dict, highsymmetrypoints)
end