using ..Paths

"""
    Lattice

Type that contains all information about a lattice.

### Fields
* `basis::Matrix`: columns are basis vectors
* `latticedim::Int`: how many of the basis vectors are lattice vectors
* `spacecoordinates::Matrix`: columns are coordinates of orbitals in the unit cell (w.r.t. basis vectors)
* `extracoordinates::Matrix`: additional non-spatial coordinates (e.g., could be sublattice index)
* `extralabels::Dict`: labels for the non-spatial coordinates
* `specialpoints::LabeledPoints`: High-symmetry points of the lattice (see `?Structure.Paths.LabeledPoints`)

### Constructors
    Lattice(basis::Matrix)
    Lattice(basis::Matrix [, latticedim::Int], spacecoordinates::Matrix [, extracoordinates::Matrix]; extralabels=Vector{String}(), specialpoints=LabeledPoints())

### Property functions
    latticedim(lat::Lattice)
    countorbitals(lat:Lattice)
    spacedim(lat::Lattice)
    extraspacedim(lat::Lattice)

    hasdimension(lat::Lattice, name::String)
    assertdimension(lat::Lattice, name::String)

    basis(lat::Lattice, ...)
    getA(lat::Lattice, ...)
    getB(lat::Lattice, ...)

    coordinates(lat::Lattice, ...)
    positions(lat::Lattice, ...)
    allpositions(lat::Lattice, ...)
    extracoordinates(lat::Lattice, ...)

    filterindices(lat::Lattice, name::String, condition::Function)

### Method functions
    setextracoordinates!
    fractionalize!
    foldfractional
    foldcoordinates!
    rotatebasis!
    rotatecoordinates!
    translate!
    displace!
    displaceZ!
    mirrorZ!
    mirrorZ
    newdimension!
    mergelattices!
    mergelattices

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
