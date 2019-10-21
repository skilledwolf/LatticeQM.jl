using SparseArrays

function graphene(lat::Lattice; V::Float64=0.0, Δ::Float64=0.0,
        λrashba=0.0, haldane=0.0, spin_orbit=0.0,
        zeeman=0.0, zeeman_staggered=0.0, zeeman_staggered_noncol=0.0,
        zeeman_layered=0.0, zeeman_layered_noncol=0.0,
        mode=:nospin, format=:auto, kwargs...)

    t(args...) = graphene_hops(args...; kwargs...)

    hops = get_hops(lat, t; format=format)

    electric_field_z!(hops, lat, V)
    sublattice_imbalance!(hops, lat, Δ)

    haldane!(hops, lat, haldane)

    if mode==:spinhalf
        hops = extend_space(hops, mode)
    end

    spin_orbit!(hops,lat,spin_orbit)
    zeeman!(hops, lat, zeeman)
    zeeman_staggered!(hops, lat, zeeman_staggered)
    zeeman_staggered_noncol!(hops, lat, zeeman_staggered_noncol)
    zeeman_layered!(hops, lat, zeeman_layered)
    zeeman_layered_noncol!(hops, lat, zeeman_layered_noncol)

    if isa(λrashba, Function) || !(λrashba≈0.0)
        H_rashba = get_bloch(rashba_hops(lat, λrashba; format=format), symmetric=false)
    else
        H_rashba = k -> 0
    end

    k -> get_bloch(hops)(k) .+ H_rashba(k)
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


cross2D(vec1, vec2) = vec1[1] * vec2[2] - vec1[2] * vec2[1]

function haldane!(hops, lat::Lattice, t2=1.0, ϕ=π/2)
    """
    This method is a somewhat inefficient way to compute the haldane hopping matrix.
    The only upside to it is that it uses methods that I already implemented and
    that it is fairly general.
    """

    if t2≈0.0
        return nothing
    end

    # NN  = find_neighbors(lat, 1.0)
    NNN = get_neighbors(lat, √3)

    N = atom_count(lat)
    R = positions(lat) # positions of atoms within unit cell
    A = get_A(lat)

    neighbors = [[i;j] for i=-1:1 for j=-1:1]
    δAs = [A * v for v in neighbors]

    # hops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()
    for δR = keys(NNN)
        if !haskey(hops, δR)
            hops[δR] = spzeros(ComplexF64, N, N)
        end
    end

    for (δR, NNN_pairs) = NNN
        for (i,j) in NNN_pairs

            Ri = R[:,i] .+ (A * δR)
            Rj = R[:,j]

            for δAi=δAs, δri=eachcol(R)

                R0 = δAi .+ δri

                if 0.9 < norm(Ri-R0) < 1.1 && 0.9 < norm(Rj-R0) < 1.1 # we found the common
                    hops[δR][i,j] += t2 * 1.0im * sign( cross2D(R0.-Rj, Ri.-R0) ) #* exp(1.0im * ϕ * sign( cross2D(R0.-Rj, Ri.-R0) ) )
                    break
                end
            end
        end
    end

    nothing
end

function spin_orbit!(hops, lat, t2, args...)

    if t2≈0.0
        return nothing
    end

    newhops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()

    haldane!(newhops, lat, t2, args...)

    newhops = kron(newhops, σZ)

    for (δR, hop)=newhops
        if !haskey(hops,δR)
            hops[δR] = hop
        else
            hops[δR] += hop
        end
    end

    nothing
end

function rashba_hops(lat::Lattice, λ::Function; format=:auto)

    function hop(r1, r2=0.0)
        δr = (r1.-r2)[1:3]
        dr = norm(δr)

        if 0.9 < dr && dr < 1.1 && abs(δr[3]) < 0.2
            δr ./= dr
            return 1.0im .* λ((r1.+r2)./2.0) .* (σX.*δr[2] .- σY.*δr[1])
        else
            return zero(σ0)
        end
    end

    get_hops(lat, hop; format=format)
end
rashba_hops(lat::Lattice, λ::Float64; kwargs...) = rashba_hops(lat, x->λ; kwargs...)

function zeeman_hops(args...; kwargs...)
    newhops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()

    zeeman!(newhops, args...; kwargs...)

    newhops
end

function zeeman!(hops, lat::Lattice, Mv::Function)
    zero0 = zeros(Int, lattice_dim(lat))

    N = atom_count(lat)
    R = positionsND(lat)

    σn(vec) = sum(vec[i] .* σs[i] for i=1:3)

    if haskey(hops, zero0)
        hops[zero0] += sum(kron(unit(i,i,N), σn(Mv(R[:,i]))) for i=1:N)
    else
        hops[zero0] = sum(kron(unit(i,i,N), σn(Mv(R[:,i]))) for i=1:N)
    end

    nothing
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

    if haskey(hops, zero0)
        hops[zero0] += kron(Matrix(1.0I, N, N), σn)
    else
        hops[zero0] = kron(Matrix(1.0I, N, N), σn)
    end

    nothing
end
zeeman!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)

function unit(i::Int,j::Int, N::Int)
    mat = zeros(ComplexF64, N, N)
    mat[i,j] = 1.0

    mat
end

function zeeman_staggered!(hops, lat::Lattice, Mv::Vector{Float64}; format=:auto)
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

    hops[zero0] += sum(kron(unit(i,i,N), sublattice[i] .* σn) for i=1:N)

    nothing
end
zeeman_staggered!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_staggered!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)


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

    hops[zero0] += sum(kron(unit(i,i,N), M0 .* σn(Mv(sublattice[i]))) for i=1:N)

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
    hops[zero0] += sum(kron(unit(i,i,N), σn(layer[i])) for i=1:N)

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

function graphene_hops(r1::T1, r2::T2; kwargs...) where {T1<:AbstractVector{Float64}, T2<:AbstractVector{Float64}}

    graphene_hops(r1.-r2; kwargs...)
end

function graphene_hops(δr::T; kwargs...) where {T<:AbstractVector{Float64}}

    graphene_hops(δr[1]^2+δr[2]^2+δr[3]^2, δr[3]^2; kwargs...)
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
