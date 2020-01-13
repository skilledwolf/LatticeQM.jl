###################################################################################################
# Properties
###################################################################################################

latticedim(lat::Lattice) = size(lat.A, 1)
countatoms(lat::Lattice) = size(lat.atoms, 2)
extraspacedim(lat::Lattice) = size(lat.atoms_aux, 1)
spacedim(lat::Lattice) = latticedim(lat) + extraspacedim(lat)

hasdimension(lat::Lattice, name::String) = haskey(lat.extradimensions, name)
assertdimension(lat::Lattice, name::String) = !hasdimension(lat, name) ? error("No $name coordinates specified.") : Nothing

getA(lat::Lattice) = lat.A
getB(lat::Lattice) = inv(transpose(lat.A))

coordinates(lat::Lattice) = lat.atoms
positions(lat::Lattice) = getA(lat) * coordinates(lat)
allpositions(lat::Lattice, args...) = [ positions(lat); extrapositions(lat, args...) ]

extrapositions(lat::Lattice) = lat.atoms_aux
extrapositions(lat::Lattice, dim::String) = extrapositions(lat, [dim])
function extrapositions(lat::Lattice, dims::Vector{String})
    for dim in aux_dim
        assertdimension(lat, dims)
    end
    indices = [lat.extradimensions[dim] for dim=dims]
    extrapositions(lat)[indices, :]
end

function setextrapositions(lat::Lattice, dim::String, values)
    assertdimension(lat, dim)
    lat.atoms_aux[lat.extradimensions[dim],:] = values
end

function filterindices(lat::Lattice, name::String, condition::Function)
    return [index for (index, val) in enumerate(extrapositions(lat, name)) if condition(val)]
end
filterpositions(lat::Lattice, args...) = positions(lat)[:,filterindices(lat, args...)]
filtercoordinates(lat::Lattice, args...) = coordinates(lat)[:,filterindices(lat, args...)]

function fractionalize(lat::Lattice, positions::AbstractMatrix{Float64})
    return inv(lat.A) * positions
end
function fractionalize!(lat::Lattice, positions::AbstractMatrix{Float64})
    positions[:] = (inv(lat.A) * positions)[:]
end
foldfractional(frac_positions::AbstractMatrix{Float64}) = mod.(frac_positions,1)

###################################################################################################
# Backwards compatibility
###################################################################################################
@legacyremoved positions3D
@legacyremoved getA_3D

export atom_count
@legacyalias countatoms atom_count

@legacyalias latticedim lattice_dim
@legacyalias extraspacedim auxspace_dim
@legacyalias spacedim space_dim

export has_dimension, assert_dimension
@legacyalias hasdimension has_dimension
@legacyalias assertdimension assert_dimension

export get_A, get_B
@legacyalias getA get_A
@legacyalias getB get_B

@legacyalias coordinates get_coordinates
@legacyalias allpositions positionsND
@legacyalias extrapositions get_positions_in
@legacyalias setextrapositions set_positions_in

@legacyalias filterindices get_filtered_indices
@legacyalias filterpositions get_filtered_positions
@legacyalias filtercoordinates get_filtered_coordinates

@legacyalias fractionalize to_fractional
@legacyalias fractionalize! to_fractional!
@legacyalias foldfractional frac_to_unitcell

