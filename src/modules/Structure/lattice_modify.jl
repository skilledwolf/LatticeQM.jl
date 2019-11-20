
function fold_atoms!(lat::Lattice)
    lat.atoms .= mod.(lat.atoms, 1)

    nothing
end

function rotate_A_XY!(lat::Lattice, α::Float64)
    @views lat.A[1:2,1:2] = rot(α) * lat.A[1:2,1:2]

    nothing
end

function rotate_atoms_XY!(lat::Lattice, θ::Float64)

    mat = inv(lat.A) * rot(θ)* lat.A

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
