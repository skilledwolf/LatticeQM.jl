import ...Utils: padvec

import LinearAlgebra: Diagonal

# Interface to Supercell module
superlattice(lat::Lattice, superperiods::Vector{Int}, args...) = superlattice(lat, Matrix(Diagonal(superperiods)), args...)
"""
    superlattice(lat::Lattice, superperiods; kwargs...)

Construct a superlattice from `lat`.

- `superperiods`: either a vector of integers (scales each direct basis vector)
  or an integer matrix `S` whose columns define the linear transformation of the
  primitive basis, i.e. the new basis is `A * S`.

Copies atoms into the new supercell and returns a new `Lattice` with updated
coordinates and `specialpoints` preserved. Extra coordinates (e.g. sublattice,
layer, z) are replicated.
"""
function superlattice(lat::Lattice, superperiods::Matrix{Int}; kwargs...)

    if superperiods==[[1,0],[0,1]] # nothing to do!
        return deepcopy(lat)
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
    superlattice(lat::Lattice, S::Matrix{Int}, supercellints::Matrix{Int}; kwargs...)

Low‑level constructor where the integer coordinates inside one supercell are
pre‑computed and provided as `supercellints` (columns). Not intended for
general users.
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


repeat(lat::Lattice, repeat=[0:0,0:0]) = repeat!(deepcopy(lat), repeat)

"""
    repeat(spacecoordinates[, Λ], repeat=[0:0, 0:0])

Tile a set of point coordinates by integer translations of the lattice. Returns
the concatenated coordinate matrix.

Arguments
- `spacecoordinates::AbstractMatrix`: columns are point positions in Cartesian
  space.
- `Λ::AbstractMatrix` (optional): direct lattice basis (defaults to 2D identity).
- `repeat::Vector{<:AbstractRange}`: ranges along a1 and a2, e.g. `[0:1, 0:1]` for 2×2.

Example
```julia
coords2x2 = repeat(coords, A, [0:1, 0:1])
```
"""
function repeat(spacecoordinates::AbstractMatrix, Λ=Matrix{Float64}(I, 2, 2), repeat=[0:0,0:0])
    @assert size(spacecoordinates,2) == 2 # only implemented for 2d lattices at the moment
    Λsuper = [Λ * [i; j] for i=repeat[1] for j=repeat[2]]
    spacecoordinates_super = hcat([spacecoordinates .+ v for v in Λsuper]...)
end

sortextraposition!(lat::Lattice, name::String) = sortposition!(lat,name) # alias for backwards compatibility
"""
    sortposition!(lat, name::String[, sortfunc])

Sort orbitals by an extra coordinate `name` (e.g. `"z"`, `"sublattice"`). Useful
for layered plots so that layers appear in front/back order. Returns the
permutation applied.
"""
function sortposition!(lat::Lattice, name::String, sortfunc=(x->x))

    perm = sortperm(sortfunc(vec(extracoordinates(lat,name))))

    lat.extracoordinates[:,:] = lat.extracoordinates[:,perm]
    lat.spacecoordinates[:,:] = lat.spacecoordinates[:,perm]

    perm
end
"""
    sortposition!(lat, index::Int[, sortfunc])

Index‑based variant of `sortposition!`, using the `index`‑th extra coordinate.
Returns the permutation applied.
"""
function sortposition!(lat::Lattice, index::Int, sortfunc=(x->x))

    perm = sortperm(sortfunc(vec(coordinates(lat,index))))

    lat.extracoordinates[:,:] = lat.extracoordinates[:,perm]
    lat.spacecoordinates[:,:] = lat.spacecoordinates[:,perm]

    perm
end


"""
    crop2unitcell!(lat)

Remove orbitals whose fractional coordinates fall outside `[0,1)` along the
primitive directions. Operates in place and returns `lat`.
"""
function crop2unitcell!(lat::Lattice)#, lat1::Lattice)
    indices = [i for (i,a) in enumerate(eachcol(lat.spacecoordinates[1:latticedim(lat),:])) if inunitrange(a)]
    lat.spacecoordinates = lat.spacecoordinates[:,indices]
    lat.extracoordinates = lat.extracoordinates[:,indices]
    lat
end
"""
    crop2unitcell(positions, Λ)

Return the subset of `positions` (Cartesian) that lie inside the primitive unit
cell spanned by `Λ` (direct basis). Returns a matrix whose columns are the kept
positions.
"""
crop2unitcell(positions::AbstractMatrix, Λ::AbstractMatrix) = crop2unitcell(inv(Λ)*positions)
function crop2unitcell(coordinates::Matrix{<:AbstractFloat})
    offset = 1e-5 * (1+0.3*rand())
    crop_iterator = filter(x->inunitrange(x; offset=offset), eachcol(coordinates))
    convert(Array, VectorOfArray(collect(crop_iterator)))
end

"""
    bistack(lat, δz; fracshift=[0.0, 0.0])

Create a bilayer by duplicating `lat` and offsetting the upper copy by `δz`
along the third (z) coordinate. Optionally shift the second layer in fractional
in‑plane coordinates by `fracshift` before stacking. Adds a `"z"` extra
coordinate if missing. Returns the new stacked `Lattice`.
"""
function bistack(lat::Lattice, δz::Float64; fracshift=[0.0; 0.0])
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

    return Lattice(A, hcat(spacecoordinates1, spacecoordinates2), hcat(extracoordinates1, extracoordinates2), lat.extralabels)
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

Given an integer lattice Z^D and an integer matrix `M` describing a (possibly
non‑orthogonal) superlattice, return the integer points inside one supercell of
the new lattice. The number of columns equals `abs(det(M))`.

The small `offset` avoids boundary ambiguities for points lying exactly on cell
faces.
"""
function supercellpoints(M::Matrix{Int}; offset::Float64=sqrt(eps()))#, check=false) # offset is a dummy variable, should be removed
    d = size(M)[1]
    n = round(Int, abs(det(M)))
    Φ = inv(M) # transformation from integer lattice into unit cube coordinates of the
    # superlattice unit cell

    # Get all corner points of the (non-orthogonal) parallepiped described by M
    # then get limits for an enclosing box region
    corners = M * UnitCubeCorners(d)
    bounds = extrema(corners; dims=2)

    innerpoints = Vector{Int}[]
    sizehint!(innerpoints, n)

    for x in Iterators.product((range(x[1]...) for x in eachrow(bounds))...)
        x = [x...]
        x0 = Φ * x
        if all((x0 .< 1 - offset) .& (x0 .>= 0 - offset))
            push!(innerpoints, x)
        end
    end

    @assert length(innerpoints)==n "Something went wrong: Could not find all cells that lie in the supercell."

    hcat(innerpoints...)
end
precompile(supercellpoints, (Matrix{Int},))



