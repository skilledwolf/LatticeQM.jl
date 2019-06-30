
################################################################################
"""
    Lattice properties
"""

lattice_dim(lat::Lattice) = size(lat.A)[1]
atom_count(lat::Lattice) = size(lat.atoms)[2]
auxspace_dim(lat::Lattice) = size(lat.atoms_aux)[1]
space_dim(lat::Lattice) = lattice_dim(lat) + auxspace_dim(lat)

has_dimension(lat::Lattice, name::String) = haskey(lat.extradimensions, name)

get_A(lat::Lattice) = lat.A

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

################################################################################
"""
    Lattice methods
"""

rot(θ::Float64) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

function fold_atoms!(lat::Lattice)
    lat.atoms .= mod.(lat.atoms, 1)

    nothing
end

function rotate_A_XY!(lat::Lattice, α::Float64)
    @views lat.A[1:2,1:2] = rot(α) * lat.A[1:2,1:2]

    nothing
end

function rotate_atoms_XY!(lat::Lattice, θ::Float64)

    mat = inv(lat.A) * [cos(θ) -sin(θ); sin(θ) cos(θ)] * lat.A

    lat.atoms[1:2,:] .= mat * lat.atoms[1:2,:]

    nothing
end

function translate!(lat::Lattice, name::String, δ::Float64)
    @assert has_dimension(lat, name)

    lat.atoms_aux[lat.extradimensions[name],:] .+= δ

    nothing
end

function mirrorZ!(lat::Lattice)
    @assert has_dimension(lat, "z")

    lat.atoms_aux[lat.extradimensions["z"],:] .*= (-1.0)

    nothing
end

function add_dimension!(lat::Lattice, name::String, aux_positions::AbstractVector{Float64})
    @assert size(lat.atoms)[2] == size(aux_positions)[1]

    lat.extradimensions[name] = auxspace_dim(lat)+1
    lat.atoms_aux = vcat(lat.atoms_aux, aux_positions')

    nothing
end

function add_dimension!(lat::Lattice, name::String, aux_positions::AbstractMatrix{Float64})
    @assert size(lat.atoms)[2] == size(aux_positions)[2]

    lat.extradimensions[name] = auxspace_dim(lat)+1
    lat.atoms_aux = vcat(lat.atoms_aux, aux_positions)

    nothing
end

function combine_lattices(lat1::Lattice, lat2::Lattice)
    """
    Note: this is the lazy option so far. A better way would be to combine at least extradimensions.
    """
    @assert lat1.A == lat2.A
    @assert lat1.extradimensions == lat2.extradimensions

    Lattice(lat1.A,
            hcat(lat1.atoms, lat2.atoms),
            hcat(lat1.atoms_aux, lat2.atoms_aux),
            lat1.extradimensions
    )
end

function combine_lattices!(lat1::Lattice, lat2::Lattice)
    """
    Note: this is the lazy option so far. A better way would be to combine at least extradimensions.
    """
    @assert lat1.A == lat2.A
    @assert lat1.extradimensions == lat2.extradimensions

    lat1.atoms = hcat(lat1.atoms, lat2.atoms)
    lat1.atoms_aux = hcat(lat1.atoms_aux, lat2.atoms_aux)

    nothing
end
