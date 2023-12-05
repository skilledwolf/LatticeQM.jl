import ..Paths

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
    addbasis!
    addorbital!
    addorbitals!
    addextra!
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
    specialpoints::Paths.LabeledPoints
end


function Lattice(basis::Matrix{Float64}, spacecoordinates::Matrix{Float64}; extralabels::Vector{String}=Vector{String}(), specialpoints=Paths.LabeledPoints())
    # @assert size(spacecoordinates,1) == size(basis,2)
    # @assert size(basis,1) == size(basis,2)

    extracoordinates = zeros(Float64, length(extralabels), size(spacecoordinates, 2)) # by default there is no perpendicular space
    Lattice(basis, spacecoordinates, extracoordinates; extralabels=extralabels, specialpoints=specialpoints)
end

import LinearAlgebra: I

function Lattice(basis::Matrix{Float64}, spacecoordinates::Matrix{Float64}, extracoordinates::Matrix{Float64}; extralabels::Vector{String}=Vector{String}(), specialpoints=Paths.LabeledPoints())
    d = size(basis, 2)
    D = size(basis, 1)
    @assert D >= d "Cannot have more basis vectors (here $d) than space dimensions (here $D)."

    if d < D
        A = Matrix(one(eltype(basis)) * I, D, D)
        A[:, 1:d] = basis[:, :]
        Lattice(A, d, spacecoordinates, extracoordinates; extralabels=extralabels, specialpoints=specialpoints)
    elseif d == D
        Lattice(basis, D, spacecoordinates, extracoordinates; extralabels=extralabels, specialpoints=specialpoints)
    end
end

import LinearAlgebra: det

function Lattice(basis::Matrix{Float64}, latticedim::Int, spacecoordinates::Matrix{Float64}, extracoordinates::Matrix{Float64}; extralabels::Vector{String}=Vector{String}(), specialpoints=Paths.LabeledPoints())
    D1 = size(basis, 1)
    D2 = size(basis, 2)
    D3 = size(spacecoordinates, 1)
    @assert D3 == D2 "Number of coordinates (here $D3) must have same length as number of basis vectors (here $D2)."
    @assert D1 == D2 "Must have same number of basis vectors (here $D2) as space dimension (here $D1)."
    @assert !isapprox(det(basis), 0; atol=sqrt(eps())) "Basis must consist of linearly-independent vectors."
    @assert D1 >= latticedim
    @assert size(extracoordinates, 2) == size(spacecoordinates, 2)

    extralabels_dict = Dict{String,UInt}(key => index for (index, key) in enumerate(extralabels))
    Lattice(basis, latticedim, spacecoordinates, extracoordinates, extralabels_dict, specialpoints)
end


function Lattice(basis::Matrix; periodic::Int=0, extra::Vector{String}=String[])
    Lattice(float(basis), periodic, zeros(Float64, size(basis, 1), 0), zeros(Float64, size(extra, 1), 0); extralabels=extra, specialpoints=Paths.LabeledPoints())
end

function Lattice()
    Lattice(zeros(0, 0))
end

function Lattice(d::Int, D::Int; kwargs...)
    Lattice(Matrix(I, D, D); periodic=D, kwargs...)
end

function addbasis!(lat::Lattice, vec::Vector, kind::Symbol=:periodic)
    @assert countorbitals(lat) == 0 "Error: at the moment basis vectors can only be added before atoms/orbitals."

    if spacedim(lat) == 0
        lat.basis = float(hcat(vec))
    else
        # @assert spacedim(lat) == size(vec,1) "Basis vector does not have the correct space dimension"

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

    lat.extralabels[label] = extraspacedim(lat) + 1

    lat
end

function addorbital!(lat::Lattice, vec::Vector)
    D = spacedim(lat) + extraspacedim(lat)
    @assert size(vec, 1) == D "Error: Coordinates must have length $D."

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


###################################################################################################
# Properties
###################################################################################################

latticedim(lat::Lattice) = lat.latticedim
countorbitals(lat::Lattice) = size(lat.spacecoordinates, 2)
extraspacedim(lat::Lattice) = length(lat.extralabels) #size(lat.extracoordinates, 1)
spacedim(lat::Lattice) = size(getA(lat), 1) # used to be latticedim(lat) + extraspacedim(lat)
allspacedim(lat::Lattice) = spacedim(lat) + extraspacedim(lat)

hasdimension(lat::Lattice, name::String) = haskey(lat.extralabels, name)
assertdimension(lat::Lattice, name::String) = !hasdimension(lat, name) ? error("No $name coordinates specified.") : nothing

basis(lat::Lattice, rselector=(:), cselector=(:)) = lat.basis[rselector, cselector]
getA(lat::Lattice, rselector=(:), cselector=(:)) = basis(lat)[:, 1:latticedim(lat)][rselector, cselector]

# (lat::Lattice)(x::AbstractMatrix) = lat.basis[:,1:size(x,1)] * x 
(lat::Lattice)(x::Vector{Int}) = (y = coordinates(lat); y[1:latticedim(lat), :] .+= x; y)
(lat::Lattice)(X::AbstractVector{Vector{Int}}) = (lat(x) for x = X)

(lat::Lattice)(x::Vector{Int}, s::Symbol) = (s == :cartesian) ? basis(lat) * lat(x) : error("Mode '$s' not known.")
(lat::Lattice)(X::AbstractVector{Vector{Int}}, s::Symbol) = (s == :cartesian) ? (basis(lat) * y for y = lat(X)) : error("Mode '$s' not known.")

"""
Calculate the dual lattice of lat.A.
Here we use the general formula \$B = A * (A^T * A)^-1\$. That works also 
when the d-dim lattice is embedded in D-dim space.
"""
function getB(lat::Lattice, rselector=(:), cselector=(:))
    d = latticedim(lat)
    A = getA(lat)[:, 1:d]
    return (A*inv(transpose(A) * A))[rselector, cselector]
end

coordinates(lat::Lattice, selector=(:)) = lat.spacecoordinates[selector, :]
allcoordinates(lat::Lattice) = [coordinates(lat); extracoordinates(lat)]

# positions(lat::Lattice) = getA(lat) * coordinates(lat)
positions(lat::Lattice, selector=(:)) = (basis(lat)*coordinates(lat))[selector, :]


allpositions(lat::Lattice, args...) = [positions(lat); extracoordinates(lat, args...)]

extracoordinates(lat::Lattice) = lat.extracoordinates
extracoordinates(lat::Lattice, dim::String) = vec(extracoordinates(lat, [dim]))
function extracoordinates(lat::Lattice, dims::Vector{String})
    for dim in dims
        assertdimension(lat, dim)
    end
    indices = [lat.extralabels[dim] for dim = dims]
    extracoordinates(lat)[indices, :]
end

function setextracoordinates!(lat::Lattice, dim::String, values)
    assertdimension(lat, dim)
    lat.extracoordinates[lat.extralabels[dim], :] = values
end

function filterindices(lat::Lattice, name::String, condition::Function)
    return [index for (index, val) in enumerate(extracoordinates(lat, name)) if condition(val)]
end
function filterindices(lat::Lattice, i::Integer, condition::Function)
    return [index for (index, val) in enumerate(coordinates(lat, i)) if condition(val)]
end
function filterindices(lat::Lattice, condition::Function)
    return [index for (index, x) in enumerate(eachcol(allpositions(lat))) if condition(x)]
end
filterpositions(lat::Lattice, args...) = positions(lat)[:, filterindices(lat, args...)]
filtercoordinates(lat::Lattice, args...) = coordinates(lat)[:, filterindices(lat, args...)]

function fractionalize(lat::Lattice, positions::AbstractMatrix{Float64})
    A = getA(lat)
    return transpose(A * inv(transpose(A) * A)) * positions
end
function fractionalize!(lat::Lattice, positions::AbstractMatrix{Float64})
    A = getA(lat)
    d = latticedim(lat)

    positions[1:d, :] .= (transpose(A * inv(transpose(A) * A)) * positions)
end
foldfractional(fracpositions::AbstractMatrix{Float64}) = mod.(fracpositions, 1)


function Base.show(io::IO, lat::Lattice)
    println(io, "Lattice dimension:     ", latticedim(lat))
    println(io, "Space dimension:       ", spacedim(lat))
    println(io, "Number of atoms:       ", countorbitals(lat))
    println(io, "Non-spatial dimension: ", extraspacedim(lat))
    println(io, "Basis:")
    show(io, basis(lat))
    println(io, "")
    println(io, "Orbital/atom coordinates: ")
    show(io, allcoordinates(lat))
end

function Base.show(io::IO, m::MIME"text/plain", lat::Lattice)
    println(io, "Lattice dimension:     ", latticedim(lat))
    println(io, "Space dimension:       ", spacedim(lat))
    println(io, "Number of atoms:       ", countorbitals(lat))
    println(io, "Non-spatial dimension: ", extraspacedim(lat))
    println(io, "Basis:")
    show(io, m, basis(lat))
    println(io, "")
    println(io, "Orbital/atom coordinates: ")
    show(io, m, allcoordinates(lat))
end


####################################################################################################
# Methods
####################################################################################################

import ..rotation2D

function uniquecoordinates!(lat::Lattice)
    lat.spacecoordinates = unique(lat.spacecoordinates; dims=2)
    lat
end

function foldcoordinates!(lat::Lattice)
    d = latticedim(lat)
    lat.spacecoordinates[1:d, :] = mod.(lat.spacecoordinates[1:d, :], 1.0)
    lat
end

function rotatebasis!(lat::Lattice, α::Float64)
    @assert spacedim(lat) >= 2
    lat.basis[1:2, :] .= (rotation2D(α) * basis(lat)[1:2, :])
    lat
end

function rotatecoordinates!(lat::Lattice, θ::Float64)
    R = positions(lat)
    D = spacedim(lat)
    d = latticedim(lat)
    R[1:2, :] .= rotation2D(θ) * R[1:2, :]
    fractionalize!(lat, R)

    lat.spacecoordinates[:, :] .= R[:, :]
    # lat.extracoordinates[1:D-d,:] .= R[d+1:D,:]
    lat
end

function translate!(lat::Lattice, name::String, δ::Float64)
    @assert hasdimension(lat, name)
    lat.extracoordinates[lat.extralabels[name], :] .+= δ
    lat
end

function translate!(lat::Lattice, i::Integer, δ::Float64)
    lat.spacecoordinates[:, :] .+= basis(lat, i) * δ
    lat
end

function translate!(lat::Lattice, v::AbstractVector)
    lat.spacecoordinates[:, :] .+= v # basis(lat) * v
    lat
end

"""
    displace!(lat, f::Function)

Function f takes orbital i at position p_i and displaces it by vector v_i = f(p_i).
"""
function displace!(lat, f::Function)

    # Modify the lattice geometry
    XY = positions(lat)

    for (i_, p) = enumerate(eachcol(XY))
        lat.spacecoordinates[:, i_] += f(p)
    end

    lat
end

"""
    displaceZ!(lat, f::Function)

Function f takes orbital i at position p_i and displaces it in the third coodinate by z_i = f(p_i).
"""
function displaceZ!(lat, f::Function)

    # Modify the lattice geometry
    XY = positions(lat)

    for (i_, p) = enumerate(eachcol(XY))
        lat.spacecoordinates[3, i_] += f(p)
    end

    lat
end

mirrorZ(lat::Lattice) = mirrorZ!(copy(lat))
function mirrorZ!(lat::Lattice)
    lat.spacecoordinates[3, :] .*= (-1.0)
    lat
end

newdimension!(lat::Lattice, name::String, extracoordinates::AbstractVector{Float64}) = newdimension!(lat, name, extracoordinates')
function newdimension!(lat::Lattice, name::String, extracoordinates::AbstractMatrix{Float64})
    @assert countorbitals(lat) == size(extracoordinates, 2)

    lat.extralabels[name] = extraspacedim(lat) + 1
    lat.extracoordinates = vcat(lat.extracoordinates, extracoordinates)
    lat
end

mergelattices(lat1::Lattice, lat2::Lattice) = mergelattices!(copy(lat1), lat2)
function mergelattices!(lat1::Lattice, lat2::Lattice)
    @assert getA(lat1) ≈ getA(lat2)
    @assert lat1.extralabels == lat2.extralabels
    lat1.spacecoordinates = hcat(lat1.spacecoordinates, lat2.spacecoordinates)
    lat1.extracoordinates = hcat(lat1.extracoordinates, lat2.extracoordinates)

    lat1
end

###############################
# Reduce dimension
###############################

function reduceto0D(lat::Lattice, N::Int)
    d = latticedim(lat)
    v = ones(Int, d)
    v[d] = N
    reduceto0D(lat, v)
end
reduceto0D(lat::Lattice, v::AbstractVector{Int}) = reduceto0D(lat, Diagonal(v))
function reduceto0D(lat::Lattice, M::AbstractMatrix{Int})
    lat1 = superlattice(lat, M)
    lat1.latticedim = 0
    kdict = Paths.LabeledPoints(
        ["Γ"],
        [[0.0]],
        ["\$\\Gamma\$"],
        ["Γ"]
    )
    lat1.specialpoints = kdict

    lat1
end


function reduceto1D(lat::Lattice, N::Int)
    d = latticedim(lat)
    v = ones(Int, d)
    v[d] = N
    reduceto1D(lat, v)
end
reduceto1D(lat::Lattice, v::AbstractVector{Int}) = reduceto1D(lat, Matrix(Diagonal(v)))
function reduceto1D(lat::Lattice, M::AbstractMatrix{Int})
    lat1 = superlattice(lat, M)
    lat1.latticedim = 1
    kdict = Paths.LabeledPoints(
        ["X1", "Γ", "X2", "Γ2"],
        [[-0.5], [0.0], [0.5], [1.0]],
        ["\$-G/2\$", "\$0\$", "\$G/2\$", "\$G\$"],
        ["X1", "Γ", "X2"]
    )
    lat1.specialpoints = kdict

    lat1
end