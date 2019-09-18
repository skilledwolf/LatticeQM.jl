
function get_operator(lat::Lattice, name::String, args...; kwargs...)

    if name == "MX" || name == "SX"
        return SX(lat, args...; kwargs...)
    elseif name == "MY" || name == "SY"
        return SY(lat, args...; kwargs...)
    elseif name == "MZ" || name == "SZ" || name == "spin"
        return SZ(lat, args...; kwargs...)
    elseif name =="Sn"
        return S_n(lat, args...; kwargs...)
    elseif name == "sublattice"
        return sublattice_N(lat, args...; kwargs...)
    elseif name == "sublatticeA"
        return sublattice_N(lat, 0, args...; kwargs...)
    elseif name == "sublatticeB"
        return sublattice_N(lat, 1, args...; kwargs...)
    elseif name == "sublatticeAspin"
        return sublattice_N(lat, 0, 2, args...; kwargs...)
    elseif name == "sublatticeBspin"
        return sublattice_N(lat, 1, 2, args...; kwargs...)
    elseif name == "Hubbard"
        return get_Hubbard(lat, args...; kwargs...)
    elseif name == "CappedYukawa"
        return get_CappedYukawa(lat, args...; kwargs...)
    else
        error("Requested operator '$name' not found.")
    end
end

get_operator(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [get_operator(lat, name, args...; kwargs...) for name=names]

get_projector(lat::Lattice, name::String, args...; kwargs...) = expval_f(get_operator(lat, name, args...; kwargs...))
get_projector(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [expval_f(get_operator(lat, name, args...; kwargs...)) for name=names]

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


################################################################################
################################################################################
################################################################################

# Scalar generators
CappedYukawa(r::AbstractVector{Float64}; kwargs...) = CappedYukawa(norm(r); kwargs...)
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
Hubbard(r::AbstractVector{Float64}; kwargs...) = Hubbard(norm(r); kwargs...)
Hubbard(r::Float64; a=0.5, U=1.0) = U * heaviside(a-r)


# Lattice operators
function get_Hubbard(lat, neighbors=[[0;0]]; mode=:nospin, format=:auto, kwargs...)
    """
    returns Dict(δL => Matrix(V(r_i-r_j+δL))_ij) where δL are vectors that
    connect unit cells. The set of δL's (in units of lattice vectors) is specified by 'neighbors'.
    """
    ee_exchange = get_hops(lat, neighbors, r->Hubbard(r; kwargs...); format=format)

    extend_space(ee_exchange, mode)
end

function get_CappedYukawa(lat, neighbors=[[i;j] for i=-1:1 for j=-1:1]; mode=:nospin, format=:auto, kwargs...)
    ee_exchange = get_hops(lat, neighbors, r->CappedYukawa(r; kwargs...); format=format)

    extend_space(ee_exchange, mode)
end

# build_CappedYukawa(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->CappedYukawa(r; kwargs...); mode=mode, format=format)
# build_Hubbard(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->Hubbard(r; kwargs...); mode=mode, format=format)
