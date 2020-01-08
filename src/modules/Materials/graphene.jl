using SparseArrays

function graphene(lat::Lattice; mode=:nospin, format=:auto, kwargs...)

    t(args...) = graphene_hops(args...; kwargs...)

    hops = get_hops(lat, t; format=format)

#     add_transversepotential!(hops, lat, V)
#     add_sublattice_imbalance!(hops, lat, Δ)

    if mode==:spinhalf
        hops = extend_space(hops, mode)
    end

#     add_spinorbit!(hops,lat,spin_orbit)
#     add_zeeman!(hops, lat, zeeman)

#     zeeman_staggered!(hops, lat, zeeman_staggered)
#     zeeman_staggered_noncol!(hops, lat, zeeman_staggered_noncol)
#     zeeman_layered!(hops, lat, zeeman_layered)
#     zeeman_layered_noncol!(hops, lat, zeeman_layered_noncol)

    hops
end


function add_sublatticeimbalance!(hops, lat::Lattice, Δ::AbstractFloat; kwargs...)

    # Only go through the trouble of constructing this matrix for finite Δ
    if abs(Δ) ≈ 0
        return nothing
    end

    μ = Δ .* (get_positions_in(lat, "sublattice") .- 0.5)
    add_chemicalpotential!(hops, lat, μ)

    nothing
end


function add_haldane!(hops, lat::Lattice, t2; ϕ=π/2, spinhalf=false, mode=:none, zmode=:none)
    """
    This method is a somewhat inefficient way to compute the haldane hopping matrix.
    The only upside to it is that it uses methods that I already implemented and
    that it is fairly general.
    """

    t2 = hcat(t2) # turn scalars into 1x1 matrix

    if norm(t2) ≈ 0
        return nothing
    end

    if spinhalf
        t2 = extend_space(t2, :spinhalf)
    end

    d = size(t2,1)

    cross2D(vec1, vec2) = vec1[1] * vec2[2] - vec1[2] * vec2[1] # needed later on in this scope

    # NN  = find_neighbors(lat, 1.0)
    NNN = get_neighbors(lat, √3)

    N = atom_count(lat)
    R = positions(lat) # positions of atoms within unit cell
    sublattice = get_positions_in(lat, "sublattice") .- 0.5
    zpositions = get_positions_in(lat, "z")
    A = get_A(lat)

    neighbors = [[i;j] for i=-1:1 for j=-1:1]
    δAs = [A * v for v in neighbors]

    # hops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()
    for δR = keys(NNN)
        if !haskey(hops, δR)
            hops[δR] = spzeros(ComplexF64, N*d, N*d)
        end
    end

    for (δR, NNN_pairs) = NNN
        for (i,j) in NNN_pairs

            Ri = R[:,i] .+ (A * δR)
            Rj = R[:,j]

            if mode==:sublatticeA
                s = (sublattice[i]>0) ? 1 : 0 #sign(sublattice[i])
            elseif mode==:sublatticeB
                s = (sublattice[i]<0) ? 1 : 0
            elseif mode==:anti
                s = sign(sublattice[i])
            else
                s = 1
            end

            if zmode==:positive
                s *= (zpositions[i]>0) ? 1 : 0 #(zpositions[i] > 0) ? 1 : 0
            elseif zmode==:negative
                s *= (zpositions[i]<0) ? 1 : 0
            elseif zmode==:anti
                s *= sign(zpositions[i])
            end

            for δAi=δAs, δri=eachcol(R)

                R0 = δAi .+ δri

                if 0.9 < norm(Ri-R0) < 1.1 && 0.9 < norm(Rj-R0) < 1.1 # we found the common
                    hops[δR][1+d*(i-1):d+d*(i-1), 1+d*(j-1):d+d*(j-1)] += s * t2 * I * 1.0im * sign( cross2D(R0.-Rj, Ri.-R0) ) #* exp(1.0im * ϕ * sign( cross2D(R0.-Rj, Ri.-R0) ) )
                    break
                end
            end
        end
    end

    nothing
end

function get_haldane_hops(args...; kwargs...)

    newhops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()

    add_haldane!(newhops, args...; kwargs...)

    newhops
end

function add_spinorbit!(hops, lat, t2, args...)

    if t2≈0.0
        return nothing
    end

    newhops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()

    add_haldane!(newhops, lat, t2, args...)

    newhops = kron(newhops, σZ)

    addhops!(hops, newhops)

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
rashba_hops(lat::Lattice, λ::AbstractFloat; kwargs...) = rashba_hops(lat, x->λ; kwargs...)

add_rashba!(hops, lat, rashba::Function; kwargs...) = addhops!(hops, rashba_hops(lat, rashba; kwargs...))
function add_rashba!(hops, lat, rashba::AbstractFloat)
    if !(rashba≈0.0)
        add_rashba!(hops, lat, x->rasbha)
    end
    nothing
end

################################################################################
################################################################################
################################################################################

"""
    This is the hopping amplitude between to twisted layers as given by overlap
    intergrals between p_z orbitals.
"""

norm2(r) = dot(r,r)
graphene_hops(r1::AbstractVector, r2::AbstractVector; kwargs...) = graphene_hops(r1.-r2; kwargs...)
graphene_hops(δr::AbstractVector; kwargs...) = graphene_hops(δr[1]^2+δr[2]^2+δr[3]^2, δr[3]^2; kwargs...)

@inline function graphene_hops(δr_sq::AbstractFloat, δz_sq::AbstractFloat;
           tz::AbstractFloat=0.46, t0::AbstractFloat=1.0, ℓinter::AbstractFloat=0.125, ℓintra::AbstractFloat=0.08,
           z::AbstractFloat=3.0, a::AbstractFloat=1.0, cutoff::AbstractFloat=5.0)
    """
        This is the hopping amplitude between to twisted layers as given by overlap
        intergrals between p_z orbitals.
        δr0_2: AbstractFloat, squared distance between points
        δz_sq: AbstractFloat, squared z-component of the distance vector
        t0: intralayer hopping amplitude
        tz: interlayper hopping amplitude t⟂
        ℓ: interlayre hopping range
        λ: intralayer hopping range
        a: intralayer atom distance
        d: interlayer distance a⟂
    """
    ℓ = ℓinter
    λ = ℓintra

    result = 0.0
    δr0 = sqrt(δr_sq)

    if !(δr0 < 1e-1 || δr0 > cutoff)

        χ = δz_sq/δr_sq

        if δz_sq < 1e-1 # intralayer hopping
            result = t0 * (1-χ) * exp(-(δr0-a)/λ)
        else # interlayer hopping
            result = tz * χ * exp(-(δr0-z)/ℓ)
        end
    end

    result
end
