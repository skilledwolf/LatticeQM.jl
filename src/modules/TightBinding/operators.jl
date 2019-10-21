
function get_operator(lat::Lattice, name::String, args...; kwargs...)

    if name == "MX" || name == "SX"
        return SX(lat, args...; kwargs...)
    elseif name == "MY" || name == "SY"
        return SY(lat, args...; kwargs...)
    elseif name == "MZ" || name == "SZ" || name == "spin"
        return SZ(lat, args...; kwargs...)
    elseif name == "spinUP"
        return Sup(lat, args...; kwargs...)
    elseif name == "spinDOWN"
        return Sdown(lat, args...;kwargs...)
    elseif name =="Sn"
        return S_n(lat, args...; kwargs...)
    elseif name =="layer"
        return layer_op(lat, args...; kwargs...)
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

get_projector(lat::Lattice, name::String, args...; kwargs...) = expvalf(get_operator(lat, name, args...; kwargs...))
get_projector(lat::Lattice, names::AbstractVector{String}, args...; kwargs...) = [expvalf(get_operator(lat, name, args...; kwargs...)) for name=names]

# Note: get_projector is a "bad" name and should be moved to get_expvalf
#       generally I use mostly for expectation values in bandstructures,
#       not projections onto basis

################################################################################
################################################################################
################################################################################

function expvalf(𝑶::AbstractMatrix)


    f(k, ψ, ϵ) = real.(ψ' * 𝑶 * ψ)

    f
end

function expvalf(𝑶::Function)

    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)

    f
end

################################################################################
################################################################################
################################################################################

function layer_op(lat::Lattice, d=1)
    """
    typically n=0 for sublattice A and n=1 for sublattice B and d=2 for spin-1/2
    """

    assert_dimension(lat, "z")

    z = get_positions_in(lat, "z")

    zmin = minimum(z)
    zmax = maximum(z)

    z = 2.0 .* ( (z .- zmin) ./ (zmax-zmin) .- 0.5 )

    # filtered_sublattice = Array{Float64}(sublattice .== n)


    kron(Diagonal(z), Diagonal(ones(d)))
end

function sublattice_N(lat::Lattice, n, d=1)
    """
    typically n=0 for sublattice A and n=1 for sublattice B and d=2 for spin-1/2
    """

    if !has_dimension(lat, "sublattice")
        error("The lattice has no sublattice defined on it.")
    end

    sublattice = get_positions_in(lat, "sublattice")

    filtered_sublattice = Array{Float64}(sublattice .== n)


    kron(Diagonal(filtered_sublattice), Diagonal(ones(d)))
end

################################################################################
################################################################################
################################################################################

function S0(lat::Lattice)
    N = atom_count(lat)

    # d.σ ⊗ 𝟙_N
    Diagonal(ones(2*N))
end

function S_n(lat::Lattice, n::Vector{Float64})
    N = atom_count(lat)

    # d.σ ⊗ 𝟙_N
    kron(Diagonal(ones(N)), sum(n[i] .* σs[i] for i=1:3))
end

SX(lat::Lattice) = S_n(lat, [1.0, 0.0, 0.0])
SY(lat::Lattice) = S_n(lat, [0.0, 1.0, 0.0])
SZ(lat::Lattice) = S_n(lat, [0.0, 0.0, 1.0])

Sup(lat::Lattice) = 0.5 .* (S0(lat) .+ SZ(lat))
Sdown(lat::Lattice) = 0.5 .* (S0(lat) .- SZ(lat))

MX = SX
MY = SY
MZ = SZ


################################################################################
################################################################################
################################################################################

macro vectorwrap(f0, N=3)
"""
This is macro is a wrapper that takes as input a function f0(x::Float64) and adds a new dispatch
f0(r1::Vector, r2::Vector) = f0(norm(r1-r2)) while making sure that r1 and r2 do not exceed length N.
"""
    return quote
        function $(esc(f0))(r1::T1, r2::Float64=0.0; kwargs...)  where {T1<:AbstractVector{Float64}}
            n = min(length(r1), $N)
            $(esc(f0))(norm(r1[1:n].-r2); kwargs...)
        end
        function $(esc(f0))(r1::T1, r2::T2; kwargs...)  where {T1<:AbstractVector{Float64}, T2<:AbstractVector{Float64}}
            n = min(length(r1), $N)
            $(esc(f0))(norm(r1[1:n].-r2[1:n]); kwargs...)
        end

        $(esc(f0))
    end
end

# Functions with scalar arguments
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))
heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
Hubbard(r::Float64; a=0.5, U=1.0) = U * heaviside(a-r)

# Functions with vector arguments
@vectorwrap CappedYukawa
@vectorwrap Hubbard

# Lattice operators
function get_Hubbard(lat, neighbors=[[0;0]]; mode=:nospin, format=:auto, kwargs...)
    """
    returns Dict(δL => Matrix(V(r_i-r_j+δL))_ij) where δL are vectors that
    connect unit cells. The set of δL's (in units of lattice vectors) is specified by 'neighbors'.
    """
    t(args...) = Hubbard(args...; kwargs...)
    ee_exchange = get_hops(lat, neighbors, t; format=format)

    extend_space(ee_exchange, mode)
end

function get_CappedYukawa(lat, neighbors=[[i;j] for i=-1:1 for j=-1:1]; mode=:nospin, format=:auto, kwargs...)
    t(args...) = Hubbard(args...; kwargs...)
    ee_exchange = get_hops(lat, neighbors, t; format=format)

    extend_space(ee_exchange, mode)
end

# build_CappedYukawa(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->CappedYukawa(r; kwargs...); mode=mode, format=format)
# build_Hubbard(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->Hubbard(r; kwargs...); mode=mode, format=format)
