import ...Utils: padvec

import LinearAlgebra: Diagonal

# Interface to Supercell module
superlattice(lat::Lattice, superperiods::Vector{Int}, args...) = superlattice(lat, Matrix(Diagonal(superperiods)), args...)
"""
    superlattice(lat::Lattice, superperiods; kwargs...)

Build superlattice geometry, i.e., build a larger supercell at it's new lattice vectors.
"""
function superlattice(lat::Lattice, superperiods::Matrix{Int}; kwargs...)

    if superperiods==[[1,0],[0,1]] # nothing to do!
        return lat
    end

    D = spacedim(lat)
    d = latticedim(lat)
    Λ = basis(lat)
    Λ[:,1:d] = Λ[:,1:d] * superperiods
    supercellints = supercellpoints(superperiods; kwargs...)

    # Existing atoms will now have to be copied
    # (note that the following code could be optimized for large systems, to use less memory -- if needed...)
    spacecoordinates_super = hcat([coordinates(lat) .+ padvec(vec,D) for vec in eachcol(supercellints)]...)
    spacecoordinates_super[1:d,:] = inv(superperiods) * spacecoordinates_super[1:d,:]
    extracoordinates_super = hcat([lat.extracoordinates for vec in eachcol(supercellints)]...) # this will copy info's such as z-coordinate or sublattice

    return Lattice(Λ, d, spacecoordinates_super, extracoordinates_super, lat.extralabels, lat.specialpoints)
end

"""
This is a low-level function. Not meant for end-users.
"""
function superlattice(lat::Lattice, superperiods::Matrix{Int}, supercellints::Matrix{Int}; kwargs...)

    D = spacedim(lat)
    d = latticedim(lat)
    Λ = basis(lat)
    Λ[:,1:d] = Λ[:,1:d] * superperiods
    # supercellints = supercellpoints(superperiods; kwargs...)

    # Existing atoms will now have to be copied
    # (note that the following code could be optimized for large systems, to use less memory -- if needed...)
    spacecoordinates_super = hcat([coordinates(lat) .+ padvec(vec,D) for vec in eachcol(supercellints)]...)
    spacecoordinates_super[1:d,:] = inv(superperiods) * spacecoordinates_super[1:d,:]
    extracoordinates_super = hcat([lat.extracoordinates for vec in eachcol(supercellints)]...) # this will copy info's such as z-coordinate or sublattice

    return Lattice(Λ, d, spacecoordinates_super, extracoordinates_super, lat.extralabels, lat.specialpoints)
end

repeat!(lat::Lattice, n::Int) = repeat!(lat, 0:n)
repeat!(lat::Lattice, ns::Vector{Int}) = repeat!(lat, [0:n for n=ns])
function repeat!(lat::Lattice, repeat::UnitRange=0:0)
    d = latticedim(lat)
    if d > 0
        repeat!(lat, collect(Iterators.repeated(repeat,d)))
    end
end

function repeat!(lat::Lattice, repeat::AbstractVector{<:AbstractRange})
    Λsuper = Vector{Int64}[]
    for I in Iterators.product(repeat...) #Iterators.repeated(1:5,3)
        append!(Λsuper, [ [I...] ])
    end

    repeat!(lat, Λsuper)
end

import ...Utils: padvec

function repeat!(lat::Lattice, intvectors::Vector{Vector{Int64}})
    lat.spacecoordinates = hcat([lat.spacecoordinates.+padvec(v,spacedim(lat)) for v in intvectors]...)
    lat.extracoordinates = hcat([lat.extracoordinates for v in intvectors]...)
    lat
end

import Base
Base.repeat(lat::Lattice, repeat=[0:0,0:0]) = repeat!(deepcopy(lat), repeat)

function repeat(spacecoordinates::AbstractMatrix, Λ=Matrix{Float64}(I, 2, 2), repeat=[0:0,0:0])
"""
Translates all points in "atom" by lattice vectors as defined by the "repeat" list of iterator.
"""
    @assert size(spacecoordinates,2) == 2 # only implemented for 2d lattices at the moment
    Λsuper = [Λ * [i; j] for i=repeat[1] for j=repeat[2]]
    spacecoordinates_super = hcat([spacecoordinates.+v for v in Λsuper]...)
end

sortextraposition!(lat::Lattice, name::String) = sortposition!(lat,name) # alias for backwards compatibility
function sortposition!(lat::Lattice, name::String, sortfunc=(x->x))
"""
Sorts all coordinates in the lattice according to extraposition "name" (e.g. z coordinates or sublattice).
This is useful for example when plotting, such that layers get plotted one on top of each other (even in supercells).
"""

    perm = sortperm(sortfunc(vec(extracoordinates(lat,name))))

    lat.extracoordinates[:,:] = lat.extracoordinates[:,perm]
    lat.spacecoordinates[:,:] = lat.spacecoordinates[:,perm]

    perm
end
function sortposition!(lat::Lattice, index::Int, sortfunc=(x->x))
"""
Sorts all coordinates in the lattice according to extraposition "name" (e.g. z coordinates or sublattice).
This is useful for example when plotting, such that layers get plotted one on top of each other (even in supercells).
"""

    perm = sortperm(sortfunc(vec(coordinates(lat,index))))

    lat.extracoordinates[:,:] = lat.extracoordinates[:,perm]
    lat.spacecoordinates[:,:] = lat.spacecoordinates[:,perm]

    perm
end


function crop2unitcell!(lat::Lattice)#, lat1::Lattice)
    indices = [i for (i,a) in enumerate(eachcol(lat.spacecoordinates[1:latticedim(lat),:])) if inunitrange(a;offset=1e-3)]
    lat.spacecoordinates = lat.spacecoordinates[:,indices]
    lat.extracoordinates = lat.extracoordinates[:,indices]
    lat
end
crop2unitcell(positions::AbstractMatrix, Λ::AbstractMatrix) = crop2unitcell(inv(Λ)*positions)
function crop2unitcell(coordinates::Matrix{<:AbstractFloat})
    offset = 1e-5 * (1+0.3*rand())
    crop_iterator = filter(x->inunitrange(x; offset=offset), eachcol(coordinates))
    convert(Array, VectorOfArray(collect(crop_iterator)))
end

function bistack(lat::Lattice, δz::Float64; fracshift=[0.0; 0.0])
    """
    Take lattice "lat", shift the copy down by δz.
    """
    if !hasdimension(lat, "z")
        N = countorbitals(lat)
        newdimension!(lat, "z", zeros(Float64,1,N))
    end

    A = getA(lat)
    spacecoordinates1 = lat.spacecoordinates
    spacecoordinates2 = deepcopy(spacecoordinates1) .+ (A*fracshift)

    extracoordinates1 = lat.extracoordinates
    extracoordinates2 = deepcopy(extracoordinates1)
    # extracoordinates1[lat.extralabels["z"],:] .+= δz/2
    extracoordinates2[lat.extralabels["z"],:] .+= δz

    return Lattice(A, hcat(spacecoordinates,spacecoordinates), hcat(extracoordinates1, extracoordinates2), lat.extralabels)
end


############################################################################################
############################################################################################
############################################################################################
############################################################################################

#### Utility functions

UnitCubeCornerIterator(dim::Int64=3) = Iterators.product(Iterators.repeated(0:1, dim)...)
UnitCubeCorners(dim::Int64=3) = hcat(([v...] for v=UnitCubeCornerIterator(dim))...)
BoundIterator(bounds) = Iterators.product(map(x->x[1]:x[2], bounds)...)
inunitrange(u; offset=sqrt(eps())) = all( (u .< 1 - offset) .& (u .>= 0 - offset) ) #all(y -> 0.0 - offset < y  < 1.0 - offset, u)

"""
    supercellpoints(M::AbstractMatrix{Int}; offset::Float64=sqrt(eps()))

Consider an integer lattice of dimension D=size(M,1). Matrix M describes a (non-orthogonal) superlattice
of this integer lattice. We want to find all lattice points that lie inside a unit cell of this
new superlattice.
"""
function supercellpoints(M::AbstractMatrix{Int}; offset::Float64=sqrt(eps()))#, check=false) # offset is a dummy variable, should be removed
    d = size(M)[1]
    Φ = inv(M) # transformation from integer lattice into unit cube coordinates of the
                # superlattice unit cell

    # Get all corner points of the (non-orthogonal) parallepiped described by M
    # then get limits for an enclosing box region
    corners = M * UnitCubeCorners(d)
    bounds = extrema(corners; dims=2)

    innerpoints = Φ * hcat(([x...] for x in BoundIterator(bounds))...)
    innerpoints = M * hcat((p for p in eachcol(innerpoints) if inunitrange(p; offset=offset))...)

    return round.(Int, innerpoints)  # M' * points
end
precompile(supercellpoints, (Matrix{Int},))

