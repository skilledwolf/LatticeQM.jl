
# Interface to Supercell module
@legacyalias superlattice build_superlattice
function superlattice(lat::Lattice, superperiods::Matrix{Int}; kwargs...)

    Λ = lat.A * superperiods
    supercellints = supercellpoints(superperiods; kwargs...)

    # Existing atoms will now have to be copied
    # (note that the following code should be optimized for large systems, to use less memory...)
    atoms_super = inv(superperiods) * hcat([lat.atoms .+ vec for vec in eachcol(supercellints)]...)
    atoms_aux_super = hcat([lat.atoms_aux for vec in eachcol(supercellints)]...) # this will copy info's such as z-coordinate or sublattice

    return Lattice(Λ, atoms_super, atoms_aux_super, lat.extradimensions, lat.highsymmetrypoints)
end

function turn0D!(lat::Lattice)
    lat.atoms = (lat.A * lat.atoms) # Turn the lattice coordinates into spatial coordinates
    lat.A = zeros(0,0) # delete the lattice vectors
    nothing
end

@legacyalias repeat! repeat_atoms!
function repeat!(lat::Lattice, repeat=[0:0,0:0])
    @assert latticedim(lat) == 2 # only implemented for 2d lattices at the moment

    Λsuper = [[i; j] for i=repeat[1] for j=repeat[2]]
    lat.atoms = hcat([lat.atoms.+v for v in Λsuper]...)
    lat.atoms_aux = hcat([lat.atoms_aux for v in Λsuper]...)
    lat
end

@legacyalias repeat repeat_atoms
repeat(lat::Lattice, repeat=[0:0,0:0]) = repeat!(deepcopy(lat), repeat)

function repeat(atoms::AbstractMatrix, Λ=Matrix{Float64}(I, 2, 2), repeat=[0:0,0:0])
"""
Translates all points in "atom" by lattice vectors as defined by the "repeat" list of iterator.
"""
    @assert size(atoms,2) == 2 # only implemented for 2d lattices at the moment
    Λsuper = [Λ * [i; j] for i=repeat[1] for j=repeat[2]]
    atoms_super = hcat([atoms.+v for v in Λsuper]...)
end


@legacyalias crop2unitcell! crop_to_unitcell!
@legacyalias crop2unitcell crop_to_unitcell
function crop2unitcell!(lat::Lattice)#, lat1::Lattice)
    indices = [i for (i,a) in enumerate(eachcol(lat.atoms)) if inunitrange(a;offset=1e-3)]
    lat.atoms = lat.atoms[:,indices]
    lat.atoms_aux = lat.atoms_aux[:,indices]
    lat
end
crop2unitcell(positions::AbstractMatrix, Λ::AbstractMatrix) = crop2unitcell(inv(Λ)*positions)
function crop2unitcell(coordinates::Matrix{<:AbstractFloat})
    offset = 1e-5 * (1+0.3*rand())
    crop_iterator = Iterators.filter(x->inunitrange(x; offset=offset), eachcol(coordinates))
    convert(Array, VectorOfArray(collect(crop_iterator)))
end

@legacyalias bistack build_bistack
function bistack(lat::Lattice, δz::Float64; fracshift=[0.0; 0.0])
    """
    Take lattice "lat", shift the copy down by δz.
    """
    if !hasdimension(lat, "z")
        N = countatoms(lat)
        newdimension!(lat, "z", zeros(Float64,1,N))
    end

    A = lat.A
    atoms1 = lat.atoms
    atoms2 = deepcopy(atoms1) .+ (A*fracshift)

    atoms_aux1 = lat.atoms_aux
    atoms_aux2 = deepcopy(atoms_aux1)
    # atoms_aux1[lat.extradimensions["z"],:] .+= δz/2
    atoms_aux2[lat.extradimensions["z"],:] .+= δz

    return Lattice(A, hcat(atoms,atoms), hcat(atoms_aux1, atoms_aux2), lat.extradimensions)
end


############################################################################################
############################################################################################
############################################################################################
############################################################################################

#### Utility functions

UnitCubeCornerIterator(dim::Int64=3) = Iterators.product(Iterators.repeated(0:1, dim)...)
UnitCubeCorners(dim::Int64=3) = convert(Array, VectorOfArray([[y...] for y=[UnitCubeCornerIterator(dim)...]]))
BoundIterator(bounds) = Iterators.product(map(x->x[1]:x[2], bounds)...)
inunitrange(u; offset=1e-6) = all( (u .< 1 - offset) .& (u .>= 0 - offset) ) #all(y -> 0.0 - offset < y  < 1.0 - offset, u)

@legacyalias supercellpoints points_within_supercell
function supercellpoints(M::AbstractMatrix{Int}; offset::Float64=1e-2)#, check=false) # offset is a dummy variable, should be removed
"""
Consider an integer lattice of dimension D=size(M)[1]. Matrix M describes a (non-orthogonal) superlattice
of this integer lattice. We want to find all lattice points that lie inside a unit cell of this
new superlattice.
"""
    d = size(M)[1]
    Φ = inv(M) # transformation from integer lattice into unit cube coordinates of the
                # superlattice unit cell
#     offset = 1e-5 + rand() * offset

    # Get all corner points of the (non-orthogonal) parallepiped described by M
    # then get limits for an enclosing box region
    corners = M * UnitCubeCorners(d)
    bounds = extrema(corners; dims=2)

    innerpoints = Φ * hcat([[x...] for x in BoundIterator(bounds)]...)
    innerpoints = M * hcat([p for p in eachcol(innerpoints) if all((p .< 1-offset) .& (p .> 0-offset))]...)

    return innerpoints  # M' * points
end


