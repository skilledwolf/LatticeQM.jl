lattice_dim(lat::Lattice) = size(lat.A)[1]
atom_count(lat::Lattice) = size(lat.atoms)[2]
auxspace_dim(lat::Lattice) = size(lat.atoms_aux)[1]
space_dim(lat::Lattice) = lattice_dim(lat) + auxspace_dim(lat)

has_dimension(lat::Lattice, name::String) = haskey(lat.extradimensions, name)

get_A(lat::Lattice) = lat.A

get_B(lat::Lattice) = inv(transpose(lat.A))

function get_A_3D(lat::Lattice)

    A = get_A(lat)
    d = lattice_dim(lat)
    N = atom_count(lat)

    diff_dim = 3 - d

    if diff_dim > 0
        # Our lattice dimension is smaller than the spatial dimension (e.g. 2D lattice in 3D space...)
        A = vcat(A, zeros(Float64, diff_dim, d) )
    end

    return A[1:3,:]
end

positions(lat::Lattice) = lat.A * lat.atoms

function positions3D(lat::Lattice)

    pos = positions(lat)
    ldim  = lattice_dim(lat)
    aux_dim = auxspace_dim(lat)
    N = atom_count(lat)

    diff_dim = 3 - ldim

    if diff_dim > 0
        # Our lattice dimension is smaller than the spatial dimension (e.g. 2D lattice in 3D space...)

        if aux_dim < diff_dim
            # We don't have enough aux-dimensions, so we'll have to pad with zeros
            pos = vcat(pos, lat.atoms_aux)
            return vcat(pos, zeros(Float64, diff_dim-aux_dim, N))
        end

        return vcat(pos, lat.atoms_aux[1:diff_dim, :])
    end

    return pos[1:3,:]
end

function get_positions_in(lat::Lattice, aux_dim::String)
    lat.atoms_aux[lat.extradimensions["sublattice"],:]
end

function get_filtered_positions(lat::Lattice, name::String, condition::Function)

    indices = [index for (index, val) in enumerate(lat.atoms_aux[lat.extradimensions[name],:]) if condition(val)]

    return positions(lat)[:,indices]
end

function get_filtered_coordinates(lat::Lattice, name::String, condition::Function)

    indices = [index for (index, val) in enumerate(lat.atoms_aux[lat.extradimensions[name],:]) if condition(val)]

    return lat.atoms[:,indices]
end

function to_fractional(lat::Lattice, positions::AbstractMatrix{Float64})
    return inv(lat.A) * positions
end
function to_fractional!(lat::Lattice, positions::AbstractMatrix{Float64})
    @views positions = inv(lat.A) * positions
end

function frac_to_unitcell(frac_positions::AbstractMatrix{Float64})
    mod.(frac_positions,1)
end
