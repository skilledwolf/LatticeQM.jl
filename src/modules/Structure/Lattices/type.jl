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

# Lattice(basis::Matrix{Float64}) = Lattice(basis, zeros(Float64, size(basis,1), 1))

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


function Lattice(basis::Matrix; periodic::Int=0, extra::Vector{String}=String[])

    Lattice(float(basis), periodic, zeros(Float64, size(basis,1), 0), zeros(Float64, size(extra,1), 0); extralabels=extra, specialpoints=kdict)
    
end

function Lattice()
    Lattice(zeros(0,0))
end

function Lattice(d::Int, D::Int; kwargs...)
    Lattice(Matrix(I, D,D); periodic=D, kwargs...)
end

function addbasis!(lat::Lattice, vec::Vector, kind::Symbol=:periodic)
    @assert countorbitals(lat) == 0 "Error: at the moment basis vectors can only be added before atoms/orbitals."

    if spacedim(lat) == 0
        lat.basis = float(hcat(vec))
    else
        @assert spacedim(lat) == size(vec,1) "Basis vector does not have the correct space dimension"
        lat.basis = hcat(lat.basis, float(vec))
    end

    if kind == :periodic
        lat.latticedim = lat.latticedim + 1
    end

    lat
end

function addextra!(lat::Lattice, label::String)
    @assert countorbitals(lat) == 0 "Error: at the moment basis can only be changed before atoms/orbitals are added."
    @assert !haskey(lat.extralabels, label) "Error: Label already exists."

    lat.extralabels[label] = extraspacedim(lat)+1

    lat
end

function addorbital!(lat::Lattice, vec::Vector)
    D = spacedim(lat)+extraspacedim(lat)
    @assert size(vec,1) == D "Error: Coordinates must have length $D."

    if countorbitals(lat) == 0
        lat.spacecoordinates = float(hcat(vec[1:spacedim(lat)]))
        lat.extracoordinates = float(hcat(vec[spacedim(lat)+1:D]))
    else
        lat.spacecoordinates = hcat(lat.spacecoordinates, float(vec[1:spacedim(lat)]))
        lat.extracoordinates = hcat(lat.extracoordinates, float(vec[spacedim(lat)+1:D]))
    end

    lat
end

function addorbitals!(lat::Lattice, orbitals::Matrix)
    for x in eachcol(orbitals)
        addorbital!(lat, Vector(x))
    end

    lat
end