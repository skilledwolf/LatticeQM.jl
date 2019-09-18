using SparseArrays

function graphene(lat::Lattice; V::Float64=0.0, Δ::Float64=0.0, zeeman::Float64=0.0, zeeman_staggered=0.0, zeeman_staggered_noncol=0.0, zeeman_layered=0.0, zeeman_layered_noncol=0.0, mode=:nospin, format=:auto, kwargs...)

    t(r::AbstractVector) = graphene_hops(r; kwargs...)

    hops = get_hops(lat, t; format=format)

    electric_field_z!(hops, lat, V)
    sublattice_imbalance!(hops, lat, Δ)

    if mode==:spinhalf

        hops = extend_space(hops, mode)

        zeeman!(hops, lat, zeeman)
        zeeman_staggered!(hops, lat, zeeman_staggered)
        zeeman_staggered_noncol!(hops, lat, zeeman_staggered_noncol)
        zeeman_layered!(hops, lat, zeeman_layered)
        zeeman_layered_noncol!(hops, lat, zeeman_layered_noncol)
    end

    get_bloch(hops)
end


function electric_field_z!(hops, lat::Lattice, V::Float64; d=3.0, kwargs...)

    # Only go through the trouble of constructing this matrix for finite V
    if abs(V) ≈ 0
        return nothing
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "z")

    N = atom_count(lat)
    R = positions(lat)

    # Get z coordinates and scale them into the unit range
    layer = get_positions_in(lat, "z")
    min = minimum(layer); max = maximum(layer)

    if abs(min-max) ≈ 0
        error("Requires at least two layers.")
    end

    layer .= (layer .- min)./(max-min)

    # get z positions and map from {0,1} to {-Δ/2,Δ/2}
    hops[zero0] += Diagonal(V .* (layer .- 0.5))

    return nothing
end


function sublattice_imbalance!(hops, lat::Lattice, Δ::Float64; kwargs...)

    # Only go through the trouble of constructing this matrix for finite Δ
    if abs(Δ) ≈ 0
        return nothing
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "sublattice")

    # diagonal matrix with on-site energies
    hops[zero0] += Diagonal(Δ .* (get_positions_in(lat, "sublattice") .- 0.5))

    return nothing
end


function zeeman!(hops, lat::Lattice, Mv::Vector{Float64}; format=:dense)
    # Only go through the trouble of constructing this matrix for finite Mv
    if norm(Mv) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "z")

    N = atom_count(lat)
    R = positions(lat)

    σn = sum(Mv[i] .* σs[i] for i=1:3)

    hops[zero0] += kron(σn, Matrix(1.0I, N, N))

    nothing
end
zeeman!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)

function unit(i::Int,j::Int, N::Int)
    mat = zeros(ComplexF64, N, N)
    mat[i,j] = 1.0

    mat
end

function zeeman_staggered(hops, lat::Lattice, Mv::Vector{Float64}; format=:auto)
    """
    Generate a Zeeman term in the Hamiltonian given a Magnetization  M(r),
    where r is a position vector.

    lat::Lattice      Lattice container
    M::Function       returns magnetization M(r)::Vector{Float} at site r
    """

    # Only go through the trouble of constructing this matrix for finite Mv
    if norm(Mv) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "sublattice")

    N = atom_count(lat)
    R = positions(lat)
    sublattice = (get_positions_in(lat, "sublattice") .- 0.5) # +Mz/2 on A and -Mz/2 on B

    σn = sum(Mv[i] .* σs[i] for i=1:3)

    hops[zero0] += sum(kron(sublattice[i] .* σn, unit(i,i,N)) for i=1:N)

    nothing
end
zeeman_staggered!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_staggered(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)


function zeeman_staggered_noncol!(hops, lat::Lattice, M0::Float64, ϕ=0.0::Float64; format=:auto)
    """
    Generate a Zeeman term in the Hamiltonian given a Magnetization  M(r),
    where r is a position vector.

    lat::Lattice      Lattice container
    M::Function       returns magnetization M(r)::Vector{Float} at site r
    """

    # Only go through the trouble of constructing this matrix for finite Mv
    if abs(M0) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "sublattice")

    N = atom_count(lat)
    R = positions(lat)
    sublattice = get_positions_in(lat, "sublattice") # +Mz/2 on A and -Mz/2 on B

    σn(Mv::T) where {T<:AbstractVector{Float64}} = sum(Mv[i] .* σs[i] for i=1:3)
    Mv(index) = [0.0,0.0,1.0] * index + [cos(ϕ), sin(ϕ), .0] * (1.0-index)

    hops[zero0] += sum(kron(M0 .* σn(Mv(sublattice[i])), unit(i,i,N)) for i=1:N)

    nothing
end


function zeeman_layered_noncol!(hops, lat::Lattice, M1::Vector{Float64}, M2::Vector{Float64})
    """
    Generate a Zeeman term in the Hamiltonian given a with a Magnetization vector
    that is
    """

    # Only go through the trouble of constructing this matrix for finite Mv
    if norm(M1) ≈ 0 && norm(M2) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "z")

    N = atom_count(lat)
    R = positions(lat)

    # Get z coordinates and scale them into the unit range
    layer = get_positions_in(lat, "z")
    min = minimum(layer)
    max = maximum(layer)

    if abs(min-max) ≈ 0
        error("Requires at least two layers.")
    end

    layer .= (layer .- min)./(max-min)

    # Define the magnetization as function of direction
    Mv(z) = M1 .* (1.0-z) + M2 .* (z)
    σn(z) = sum(Mv(z)[i] .* σs[i] for i=1:3)

    # hops[zero0] += blockdiag([sparse(σn(z_scaled)) for z_scaled=layer]...)
    hops[zero0] += sum(kron(σn(layer[i]), unit(i,i,N)) for i=1:N)

    nothing
end
zeeman_layered_noncol!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_layered_noncol!(hops, lat::Lattice, [M0,0.0,0.0], [0.0,0.0,M0]; kwargs...)
zeeman_layered!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_layered_noncol!(hops, lat::Lattice, [0.0,0.0,M0], [0.0,0.0,-M0]; kwargs...)#zeeman_field_layered!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)


################################################################################
################################################################################
################################################################################

"""
    This is the hopping amplitude between to twisted layers as given by overlap
    intergrals between p_z orbitals.
"""

norm2(r) = dot(r,r)

function graphene_hops(r1::T, r2::T; kwargs...) where {T<:AbstractVector{Float64}}

    graphene_hops(r1.-r2; kwargs...)
end

function graphene_hops(δr::T; kwargs...) where {T<:AbstractVector{Float64}}

    graphene_hops(norm2(δr), δr[3]^2; kwargs...)
end

@inline function graphene_hops(δr_sq::Float64, δz_sq::Float64;
           tm::Float64=0.46, t0::Float64=1.0, ℓ::Float64=0.125, λ::Float64=0.08,
           z::Float64=3.0, a::Float64=1.0, cutoff::Float64=5.0)
    """
        This is the hopping amplitude between to twisted layers as given by overlap
        intergrals between p_z orbitals.
        δr0_2: Float64, squared distance between points
        δz_sq: Float64, squared z-component of the distance vector
        t0: intralayer hopping amplitude
        tm: interlayper hopping amplitude t⟂
        ℓ: interlayre hopping range
        λ: intralayer hopping range
        a: intralayer atom distance
        d: interlayer distance a⟂
    """

    result = 0.0
    δr0 = sqrt(δr_sq)

    if !(δr0 < 1e-1 || δr0 > cutoff)

        χ = δz_sq/δr_sq

        if δz_sq < 1e-1 # intralayer hopping
            result = t0 * (1-χ) * exp(-(δr0-a)/λ)
        else # interlayer hopping
            result = tm * χ * exp(-(δr0-z)/ℓ)
        end
    end

    result
end
