using SparseArrays

function graphene(lat::Lattice; mode=:nospin, format=:auto, kwargs...)

    t(args...) = t_graphene(args...; kwargs...)

    hops = gethops(lat, t; format=format)

    if mode==:spinhalf
        hops = addspin(hops, mode)
    end

    hops
end

valleyoperator(args...; kwargs...) = valleyoperator!(Hops(), args...; kwargs...)
function valleyoperator!(hops, lat::Lattice; spinhalf=false, zmode=:anti, kwargs...)
    @assert latticedim(lat) == 2 && countorbitals(lat) > 1
    addhaldane!(hops, lat, √3/9; spinhalf=spinhalf, mode=:anti, zmode=zmode, kwargs...)
end

@legacyalias addsublatticeimbalance! add_sublatticeimbalance
function addsublatticeimbalance!(hops, lat::Lattice, Δ::AbstractFloat; kwargs...)

    # Only go through the trouble of constructing this matrix for finite Δ
    if abs(Δ) ≈ 0
        return nothing
    end

    μ = Δ .* (extrapositions(lat, "sublattice") .- 0.5)
    addchemicalpotential!(hops, lat, vec(μ))

    nothing
end

@legacyalias addhaldane add_haldane
function addhaldane!(hops, lat::Lattice, t2; ϕ=π/2, spinhalf=false, mode=:none, zmode=:none)
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
        t2 = addspin(t2, :spinhalf)
    end

    d = size(t2,1)

    cross2D(vec1, vec2) = vec1[1] * vec2[2] - vec1[2] * vec2[1] # needed later on in this scope

    # NN  = find_neighbors(lat, 1.0)
    NNN = Structure.getneighbors(lat, √3)

    N = countorbitals(lat)
    R = positions(lat) # positions of atoms within unit cell
    sublattice = extrapositions(lat, "sublattice") .- 0.5
    if hasdimension(lat, "z")
        zpositions = extrapositions(lat, "z")
    end
    A = getA(lat)

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

            if hasdimension(lat, "z")
                if zmode==:positive
                    s *= (zpositions[i]>0) ? 1 : 0 #(zpositions[i] > 0) ? 1 : 0
                elseif zmode==:negative
                    s *= (zpositions[i]<0) ? 1 : 0
                elseif zmode==:anti
                    s *= sign(zpositions[i]+1e-4)
                end
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

    hops
end

@legacyalias gethaldane get_haldane_hops
function gethaldane(args...; kwargs...)
    newhops = Hops()
    addhaldane!(newhops, args...; kwargs...)
end

@legacyalias addspinorbit! add_spinorbit!
function addspinorbit!(hops, lat, t2, args...)

    if t2≈0.0
        return nothing
    end

    newhops = Hops()
    addhaldane!(newhops, lat, t2, args...)
    addhops!(hops, kron(newhops, σZ))
end

@legacyalias getrashba rashba_hops
function getrashba(lat::Lattice, λ::Function; format=:auto)

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

    gethops(lat, hop; format=format)
end
getrashba(lat::Lattice, λ::AbstractFloat; kwargs...) = getrashba(lat, x->λ; kwargs...)

@legacyalias addrashba! add_rashba!
addrashba!(hops, lat, rashba::Function; kwargs...) = addhops!(hops, getrashba(lat, rashba; kwargs...))
function addrashba!(hops, lat, rashba::AbstractFloat)
    if !(rashba≈0.0)
        addrashba!(hops, lat, x->rasbha)
    end
end

################################################################################
################################################################################
################################################################################

"""
    This is the hopping amplitude between to twisted layers as given by overlap
    intergrals between p_z orbitals.
"""

norm2(r) = dot(r,r)
t_graphene(r1::AbstractVector, r2::AbstractVector; kwargs...) = t_graphene(r1.-r2; kwargs...)
t_graphene(δr::AbstractVector; kwargs...) = t_graphene(δr[1]^2+δr[2]^2+δr[3]^2, δr[3]^2; kwargs...)

@legacyalias t_graphene graphene_hops
@inline function t_graphene(δr_sq::AbstractFloat, δz_sq::AbstractFloat;
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
