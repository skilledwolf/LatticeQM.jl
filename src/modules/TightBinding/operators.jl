
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
    mat = spzeros(Complex, 2*N, 2*N)
    σn = sum(n[i] .* σs[i] for i=1:3)

    @simd for i = 1:2:2*N
        mat[i:i+1, i:i+1] .= σn
    end

    mat
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
heaviside(x) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
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


###################################################################################################
###################################################################################################
###################################################################################################

function set_filling!(hops, lat, filling; nk=100, kwargs...)
    kgrid = regulargrid(nk=nk)
    set_filling!(hops, lat, kgrid, filling; kwargs...)
end

function set_filling!(hops, lat, kgrid, filling; T=0.0)
    hops0 = get_dense(hops)
    μ = chemical_potential(get_bloch(hops0), kgrid, filling; T=T)
    add_chemicalpotential!(hops, lat, -μ)

    μ
end

function add_chemicalpotential!(hops, lat::Lattice, μ::T; localdim::Int=-1) where T<:AbstractVector{<:Float64}
    zero0 = zeros(Int, lattice_dim(lat))
    N = atom_count(lat)

    if localdim < 0 # if localdim is not set, we determine it from matrix dimensions
        D = hopdim(hops)
        d = div(D, N)
    else
        d = localdim
    end

    @assert N == length(μ)

    newhops = Dict( zero0 => sparse( (1.0+0.0im).* Diagonal(kron(μ, ones(d))) ) )
    add_hoppings!(hops, newhops)

    nothing
end

function add_chemicalpotential!(hops, lat::Lattice, μ::Function; kwargs...)
    R = positionsND(lat)
    add_chemicalpotential!(hops, lat, [μ(r) for r=eachcol(R)]; kwargs...)

    nothing
end
add_chemicalpotential!(hops, lat::Lattice, μ::Float64; kwargs...) = add_chemicalpotential!(hops, lat, μ.*ones(atom_count(lat)); kwargs...)


###################################################################################################
###################################################################################################
###################################################################################################


function get_mf_functional_new(hops, v::Dict{Vector{Int},T2}) where {T1<:Complex, T2<:AbstractMatrix{T1}}
    """
        This method takes the Hamiltonian single-particle operator h and an
        interaction potential v and returns mean-field functionals
            ℋ, E  s.t.  h_mf = ℋ[ρ]  and  ϵ_scalar = E[ρ].

        These functionals can be used to search for a self-consistent solution
        using solve_selfconsistent(...).
    """

    mf_hops, E = get_mf_hops(v)
    H(ρ) = get_bloch(get_dense(add_hoppings(hops, mf_hops(ρ))))

    H, E
end

function get_mf_hops(v::Dict{Vector{Int},T2}) where {T1<:Complex, T2<:AbstractMatrix{T1}}
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ,k] that builds the mean field hamiltonian
        (i.e. h_v(k) = 𝒱[ρ,k]).

        This may look harmless but requires a careful derivation.
    """

    # d = size(first(values(v)),1)
    # vsym(L::Vector{Int}) = 0.5 .* (v[L].+(v[L])')
    V0 = sum(v[L] for L in keys(v))

    function mf_op(ρs::Dict{Vector{Int},T2}) where {T1<:Complex, T2<:AbstractMatrix{T1}}
        δL0 = zero(first(keys(ρs)))

        # fock
        hops_mf = Dict(δL => v[δL] .* ρL for (δL,ρL) in ρs)
        # hartree
        hops_mf[δL0] .+= spdiagm(0 => V0 * diag(ρs[δL0]))

        hops_mf
    end

    function mf_scalar(ρs::Dict{Vector{Int},T2}) where {T1<:Complex, T2<:AbstractMatrix{T1}}
        δL0 = zero(first(keys(ρs)))

        # Hartree contribution
        vρ = diag(ρs[δL0])
        e_hartree = - 1/2 * (vρ' * V0 * vρ)
        @assert imag(e_hartree) ≈ 0

        # Fock contribution
        e_fock =  1/2 * sum(sum(ρL .* conj.(ρL) .* v[L] for (L,ρL) in ρs))
        @assert imag(e_hartree) ≈ 0

        real(e_hartree + e_fock)
    end

    mf_op, mf_scalar
end

###################################################################################################
###################################################################################################
###################################################################################################

function initial_guess(v::Dict{Vector{Int},T2}, mode=:random; lat=:nothing) where {T1<:Complex, T2<:AbstractMatrix{T1}}
    N = size(first(values(v)), 1)

    ρs = Dict{Vector{Int},Matrix{ComplexF64}}()
    for δL=keys(v)
        ρs[δL] = zeros(ComplexF64, size(v[δL]))
    end

    if mode==:randombig
        mat = rand(ComplexF64, N, N)
        ρs[zero(first(keys(ρs)))] = (mat + mat') ./ 2

    elseif mode==:random
        @assert mod(N,2)==0
        n = div(N,2)
        mat = rand(ComplexF64, n, n)

        # Generate a random spin orientation at a lattice site
        function randmat()
            # d = rand(Float64, 3)
            d = -1.0 .+ 2 .* rand(Float64, 3)
            p = 0.5 .* (σ0 .+ sum(d[i_]/norm(d) .* σs[i_] for i_=1:3))
        end

        ρs[zero(first(keys(ρs)))] = Matrix(sum(kron(sparse([i_],[i_], [1.0+0.0im], n,n), randmat()) for i_=1:n))

    elseif mode==:antiferro || mode==:antiferroZ
        sublA, sublB = get_operator(lat, ["sublatticeA", "sublatticeB"])
        mat = sublA .- sublB

        σUP = 0.5 .* (σ0 .+ σZ)

        ρs[zero(first(keys(ρs)))] = kron(mat,σUP)

    elseif mode==:antiferroX
        sublA, sublB = get_operator(lat, ["sublatticeA", "sublatticeB"])
        mat = sublA .- sublB

        σUP = 0.5 .* (σ0 .+ σX)

        ρs[zero(first(keys(ρs)))] = kron(mat,σUP)

    elseif mode==:ferro || mode==:ferroZ #|| mode==:ferroz
        @assert mod(N,2)==0
        n = div(N,2)
        σUP = 0.5 .* (σ0 .+ σZ)
        ρs[zero(first(keys(ρs)))] =  2. * kron(Diagonal(ones(n)),σUP)

    elseif mode==:ferroX #|| mode==:ferroz
        @assert mod(N,2)==0
        n = div(N,2)
        σLEFT = 0.5 .* (σ0 .+ σX)
        ρs[zero(first(keys(ρs)))] =  2. * kron(Diagonal(ones(n)), σLEFT)
    # elseif mode==:ferrox
    #     @assert mod(N,2)==0
    #     n = div(N,2)
    #     ρs[zero(first(keys(ρs)))] = kron(σX, Diagonal(ones(n)))

    else
        error("Unrecognized mode '$mode' in initialize_ρ(...).")

    end

    ρs
end