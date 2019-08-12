
function get_operator(lat::Lattice, name::String)

    if name == "MX" || name == "SX"
        return SX(lat)
    elseif name == "MY" || name == "SY"
        return SY(lat)
    elseif name == "MZ" || name == "SZ" || name == "spin"
        return SZ(lat)
    elseif name == "sublatticeA"
        return sublattice_N(lat, 0)
    elseif name == "sublatticeB"
        return sublattice_N(lat, 1)
    elseif name == "sublatticeAspin"
        return sublattice_N(lat, 0, 2)
    elseif name == "sublatticeBspin"
        return sublattice_N(lat, 1, 2)
    else
        error("Requested operator '$name' not found.")
    end
end

get_operator(lat::Lattice, names::AbstractVector{String}) = [get_operator(lat, name) for name=names]

get_projector(lat::Lattice, name::String) = expval_f(get_operator(lat, name))
get_projector(lat::Lattice, names::AbstractVector{String}) = [expval_f(get_operator(lat, name)) for name=names]

################################################################################
################################################################################
################################################################################

function expval_f(𝑶::AbstractMatrix)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶 * ψ)

    f
end

################################################################################
################################################################################
################################################################################

function sublattice_N(lat::Lattice, n, d=1)
    """
    typically n=0 for sublattice A and n=1 for sublattice B and d=2 for spin-1/2
    """

    if !has_dimension(lat, "sublattice")
        error("The lattice has no sublattice defined on it.")
    end

    sublattice = get_positions_in(lat, "sublattice")

    filtered_sublattice = Array{Float64}(sublattice .== n)


    kron(Diagonal(ones(d)), Diagonal(filtered_sublattice))
end

################################################################################
################################################################################
################################################################################

function S_n(lat::Lattice, n::Vector{Float64})
    N = atom_count(lat)

    kron( sum(n[i] .* σs[i] for i=1:3) , Diagonal(ones(N)))
end

SX(lat::Lattice) = S_n(lat, [1.0, 0.0, 0.0])
SY(lat::Lattice) = S_n(lat, [0.0, 1.0, 0.0])
SZ(lat::Lattice) = S_n(lat, [0.0, 0.0, 1.0])

MX = SX
MY = SY
MZ = SZ
