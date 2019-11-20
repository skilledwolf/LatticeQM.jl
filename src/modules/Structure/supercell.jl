
# Interface to Supercell module
function build_superlattice(lat::Lattice, superperiods::Matrix{Int}; kwargs...)

    Λ = lat.A * superperiods
    superlattice_int = points_within_supercell(superperiods; kwargs...)

    # Existing atoms will now have to be copied
    # (note that the following code should be optimized for large systems, to use less memory...)
    atoms_super = inv(superperiods) * hcat([lat.atoms .+ vec for vec in eachcol(superlattice_int)]...)
    atoms_aux_super = hcat([lat.atoms_aux for vec in eachcol(superlattice_int)]...) # this will copy info's such as z-coordinate or sublattice

    return Lattice(Λ, atoms_super, atoms_aux_super, lat.extradimensions, lat.highsymmetrypoints)
end

function turn0D!(lat::Lattice)
    lat.atoms = (lat.A * lat.atoms) # Turn the lattice coordinates into spatial coordinates
    lat.A = zeros(0,0) # delete the lattice vectors
    nothing
end

function repeat_atoms!(lat::Lattice, repeat=[0:0,0:0])

    Λsuper = [[i; j] for i=repeat[1] for j=repeat[2]]

    lat.atoms = hcat([lat.atoms.+v for v in Λsuper]...)
    lat.atoms_aux = hcat([lat.atoms_aux for v in Λsuper]...)

    nothing
end
function repeat_atoms(lat::Lattice, repeat=[0:0,0:0])
    lat2 = deepcopy(lat)
    repeat_atoms!(lat2, repeat)
    return lat2
end

function crop_to_unitcell!(lat::Lattice)#, lat1::Lattice)
    atoms = lat.atoms

    indices = [i for (i,a) in enumerate(eachcol(atoms)) if inunitrange(a;offset=1e-3)]

    lat.atoms = lat.atoms[:,indices]
    lat.atoms_aux = lat.atoms_aux[:,indices]

end

function build_bistack(lat::Lattice, δz::Float64; fracshift=[0.0; 0.0])
    """
    Take lattice "lat", shift the copy down by δz.
    """
    if !has_dimension(lat, "z")
        N = atom_count(lat)
        add_dimension!(lat, "z", zeros(Float64,1,N))
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


function points_within_supercell(M::AbstractMatrix{Int}; offset::Float64=1e-2)#, check=false) # offset is a dummy variable, should be removed
"""
Consider an integer lattice of dimension D=size(M)[1]. M describes a (non-orthogonal) superlattice
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

#### OLD BUGGY VERSION!!
# function points_within_supercell(M::AbstractMatrix{Int}; offset::Float64=1e-2)#, check=false)
# """
# Consider an integer lattice of dimension D=size(M)[1]. M describes a (non-orthogonal) superlattice
# of this integer lattice. We want to find all lattice points that lie inside a unit cell of this
# new superlattice.
# """
#     d = size(M)[1]
#     Φ = inv(M) # transformation from integer lattice into unit cube coordinates of the
#                 # superlattice unit cell
#     offset = 1e-5 + rand() * offset
#
#     # Get all corner points of the (non-orthogonal) parallepiped described by M
#     # then get limits for an enclosing box region
#     corners = M * UnitCubeCorners(d)
#     bounds = extrema(corners; dims=2)
#
#     LatticeIterator = Base.Generator(x -> [x...], BoundIterator(bounds)) # Iterate over all points in that lie in the enclosing box region
#     PointsIterator = Iterators.filter(x->inunitrange(Φ * x, offset=offset), LatticeIterator) # Iterator that filters out elements that lie inside the parallepiped of M
#
#     # Pre-allocate memory
#     points = Array{Int}(undef, d, round(Int, abs(det(M))))
#     for (i,point) in enumerate(PointsIterator)
#         points[:,i] .= point#round.(Int, M * point)
#     end
#
# #         points = convert(Array, VectorOfArray(collect(PointsIterator))) # Evaluate the iterators
# #         if check & (round(Int, abs(det(M))) != size(points)[2]) # Check if we found alle points
# #             throw("Oops. Did not find all lattice points of the supercell. Possibly an implementation error.")
# #         end
#
#     return points  # M' * points
# end

precompile(points_within_supercell,(AbstractMatrix{Int},))

function repeat_atoms(atoms::Matrix{<:AbstractFloat}, Λ=Matrix{Float64}(I, 2, 2), repeat=[0:0,0:0])
"""
Translates all points in "atom" by lattice vectors as defined by the "repeat" list of iterator.
atoms: Matrix with columns as point coordinates
Λ: Matrix with columns as lattice vectors
repeat: listing of iterators for each lattice vector in Λ

return:  matrix with columns as point coordinates

This is procedure is often used when constructing supercells.
"""
#         Λ = [Λ zeros(2); zeros(1,2) 1] # extend to 3D
    Λsuper = [Λ * [i; j] for i=repeat[1] for j=repeat[2]]

    atoms_super = hcat([atoms.+v for v in Λsuper]...)
end

function crop_to_unitcell(atoms::Matrix{<:AbstractFloat}, Λ=Matrix{Float64}(I, 2, 2))
"""
Returns all points that lie inside the unit cell.

atoms: Matrix with columns as point coordinates
Λ: Matrix with columns as lattice vectors

return:  matrix with columns as point coordinates

This procedure can be useful when reshaping unit cells and supercells.
"""

    Φ = inv(Λ)

    randnum = rand()
    offset = 1e-3 * (1+0.3*randnum)

    crop_iterator = Iterators.filter(x->inunitrange(Φ*x; offset=offset), eachcol(atoms))

    convert(Array, VectorOfArray(collect(crop_iterator)))
end
