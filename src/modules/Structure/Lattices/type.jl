using ..Paths

"""
    Lattice

Type that contains all information about a lattice, such as lattice vectors
in ``Lattice.a``, fractional coordinates in ``Lattice.coordinates`` and
extra dimensions in ``Lattice.extralabels``.

# Fields
- `a::Matrix{Float64}` :
  :math:`d \times d` matrix of lattice vectors a[:,i]
- `coordinates::Matrix{Float64}` :
  :math:`D \times d` matrix of fractional coordinates x[1:d,i]
  The coordinates x[d+1:D, i] are additional coordinates
- `extralabels::Dict{String,Int}` :
  Keys are the names of the additional coordinates, e.g. "sublattice" or "z".
  The values are indices for the corresponding row in ``Lattice.coordinates``.
"""
mutable struct Lattice
    basis::Matrix{Float64} # D × D
    latticedim::Int
    spacecoordinates::Matrix{Float64} # d × N # fractional coordinates w.r.t. A
    extracoordinates::Matrix{Float64} # D × N
    extralabels::Dict{String,Int}
    specialpoints::LabeledPoints
end

kdict = LabeledPoints(Dict{String,AbstractVector{Float64}}(), Dict{String,String}(), Vector{String}([]))

Lattice(basis::Matrix{Float64}) = Lattice(basis, zeros(Float64, size(basis,1), 1))

function Lattice(basis::Matrix{Float64}, spacecoordinates::Matrix{Float64}; extralabels::Vector{String}=Vector{String}(), specialpoints=kdict)
    # @assert size(spacecoordinates,1) == size(basis,2)
    # @assert size(basis,1) == size(basis,2)

    extracoordinates = zeros(Float64, length(extralabels), size(spacecoordinates,2)) # by default there is no perpendicular space
    Lattice(basis, spacecoordinates, extracoordinates; extralabels=extralabels, specialpoints=specialpoints)
end

function Lattice(basis::Matrix{Float64}, spacecoordinates::Matrix{Float64}, extracoordinates::Matrix{Float64}; extralabels::Vector{String}=Vector{String}(), specialpoints=kdict)
    d = size(basis,2)
    D = size(basis,1)
    @assert D >= d "Cannot have more basis vectors (here $d) than space dimensions (here $D)."

    if d < D
      A = Matrix(one(eltype(basis))*I, D, D)
      A[:, 1:d] = basis[:,:]  
      Lattice(A, d, spacecoordinates, extracoordinates; extralabels=extralabels, specialpoints=specialpoints)
    elseif d == D
      Lattice(basis, D, spacecoordinates, extracoordinates; extralabels=extralabels, specialpoints=specialpoints)
    end
end

function Lattice(basis::Matrix{Float64}, latticedim::Int, spacecoordinates::Matrix{Float64}, extracoordinates::Matrix{Float64}; extralabels::Vector{String}=Vector{String}(), specialpoints=kdict)
    D1 = size(basis,1); D2 = size(basis,2)
    D3 = size(spacecoordinates,1)
    @assert D3 == D2 "Number of coordinates (here $D3) must have same length as number of basis vectors (here $D2)."
    @assert D1 == D2 "Must have same number of basis vectors (here $D2) as space dimension (here $D1)."
    @assert !isapprox(det(basis), 0; atol=sqrt(eps())) "Basis must consist of linearly-independent vectors."
    @assert D1 >= latticedim
    @assert size(extracoordinates,2) == size(spacecoordinates,2)

    extralabels_dict = Dict{String,UInt}(key=>index for (index,key) in enumerate(extralabels))
    Lattice(basis, latticedim, spacecoordinates, extracoordinates, extralabels_dict, specialpoints)
end
