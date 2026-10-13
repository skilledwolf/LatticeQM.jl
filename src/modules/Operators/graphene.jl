using SparseArrays

import ..TightBinding: Hops

"""
    graphene(lat::Lattice; vectorized=true, mode=:nospin, format=:auto, cellrange=2, kwargs...)

Construct a graphene tight‑binding Hamiltonian on lattice `lat`.

- `mode`: `:nospin` for spinless, `:spinhalf` to attach a spin‑1/2 degree of freedom.
- `format`: storage for the resulting hopping object (`:auto`, `:dense`, or `:sparse`).
- `cellrange`: neighbour shell used when populating hoppings.

Returns a `Hops` object that can be passed to spectrum routines.
"""

function graphene(lat::Lattice; vectorized=true, mode=:nospin, format=:auto, cellrange=2, kwargs...)

    t(args...) = t_graphene(args...; kwargs...)

    hops = Hops(lat, t; cellrange=cellrange, format=format, vectorized=vectorized)

    if mode==:spinhalf
        hops = addspin(hops, mode)
    end

    hops
end

function addgraphene!(hops, lat::Lattice; kwargs...)
    addhops!(hops, graphene(lat; kwargs...))
end

# default parameters taken from PHYSICAL REVIEW B 82, 035409 (2010), Table 1
function graphene_rhombohedral(lat; spin=false, a=1.0, d=3.0, γ0=-3.16, γ1=0.502, γ2=-0.0171, γ3=-0.377, γ4=-0.099, kwargs...)

    function t(r1, r2=0.0)
        δr=r1.-r2
        
        if abs(norm(δr[1:2])-a)<0.01 && abs(δr[3]) < 0.01 # same layer, NN
            return γ0
        elseif abs(norm(δr[1:2]))<0.01 && abs(abs(δr[3])-d) < 0.01  ## NN-layer, vertical (dimer)
            return γ1
        elseif abs(norm(δr[1:2]))<0.01 && abs(abs(δr[3])-2*d) < 0.01 ## NNN-layer, vertical
            return γ2
        elseif abs(norm(δr[1:2])-a)<0.01 && abs(abs(δr[3])-d) < 0.01 && abs(δr[4]) < 0.01 ## NN-layer, same sublattice (A1-A2 / B1-B2)
            return γ4
        elseif abs(norm(δr[1:2])-a)<0.01 && abs(abs(δr[3])-d) < 0.01 && abs(δr[4]) > 0.01 ## NN-layer, non-dimer A1-B2: in-plane distance is a, like γ4 (the old 2a selector tagged a spurious far shell and missed every physical γ3 bond)
            return γ3
        end
        
        return 0.0
    end

    hops = Hops(lat, t; vectorized=false, kwargs...)

    if spin
        hops = addspin(hops, :spinhalf)
    end

    hops
end

import ..Structure.Lattices

"""
    addsublatticeimbalance!(hops, lat, Δ; kwargs...)

Add a sublattice‑staggered chemical potential (imbalance) of magnitude `Δ` to
`hops` on lattice `lat`. The on‑site energy is shifted by
`Δ*(s - 1/2)` with sublattice coordinate `s` (A: `s=0`, B: `s=1`), so positive
`Δ` lowers sublattice A by `Δ/2` and raises sublattice B by `Δ/2`.

No‑op for `Δ≈0`.
"""
function addsublatticeimbalance!(hops, lat::Lattice, Δ::Real; kwargs...)

    # Only go through the trouble of constructing this matrix for finite Δ
    if isapprox(Δ,0; atol=sqrt(eps()))
        return nothing
    end

    μ = Δ .* (Lattices.extracoordinates(lat, "sublattice") .- 0.5)
    addchemicalpotential!(hops, lat, vec(μ))

    nothing
end

"""
    valley(lat; spinhalf=false, kwargs...)

Construct a valley operator on `lat`. If `spinhalf=true`, attach spin structure
to match spinful Hamiltonians.
"""
function valley(lat, args...; spinhalf=false, kwargs...) 
    h=Hops()
    addvalley!(h, lat, args...; kwargs...)

    if spinhalf
        h = addspin(h, :spinhalf)
    end
    
    autoconversion(h, Lattices.countorbitals(lat))
end

"""
    addvalley!(hops, lat, fz=x->sign(x[3]+1e-3); kwargs...)

Add a valley mass term to `hops` on `lat`. Customise the layer sign via `fz`.
"""
function addvalley!(hops, lat::Lattice, fz::Function=x->sign(x[3]+1e-3); kwargs...)
    @assert Lattices.countorbitals(lat) > 1 #latticedim(lat) == 2

    t20 = √3/9
    f(R) = sign(R[4]-0.5)
    addhaldane!(hops, lat, x-> t20 * f(x) * fz(x); kwargs...)
end

include("haldane.jl")

function addspinorbit!(hops, lat, t2, args...)

    if isapprox(t2, 0; atol=sqrt(eps()))
        return nothing
    end

    newhops = Hops()
    addhaldane!(newhops, lat, t2, args...)
    addhops!(hops, kron(newhops, σZ))
end

using ..TightBinding: MAX_DENSE, MAX_DIAGS

function getrashba(lat::Lattice, λ::Function; cellrange=2, format=:auto, tmin=1e-7)

    function hop(R1, R2)
        # this is a "vectorized" hopping function, i.e., it expects
        # coordinate matrices R1 and R2 as inputs
        # (had to be implemented this way for performance)

        N = size(R1,2)
        d = 2

        # Preallocate memory: important for huge sparse matrices
        maxind = (N>MAX_DENSE) ? round(Int, MAX_DIAGS * N) : N^2 # MIN_SPARSITY * N^2 # semi-arbitrary limit for dense allocation
        IS = Vector{Int}(undef, maxind*d^2)
        JS = similar(IS)
        VS = similar(IS, ComplexF64)
        δR = similar(R1)
        R0 = similar(R1)

        count = 0
        # do not tamper with this loop.
        # this way is orders of magnitude faster than any other implementation
        @fastmath @inbounds for j=1:N
            @views δR .= R1 .- R2[:,j]
            @views R0 .= (R1 .+ R2[:,j])./2

            for i=1:N
                @views Δ = sqrt(sum(abs2,δR[1:3,i]))

                if 1.1 < Δ || Δ < 0.9 || abs(δR[3,i]) > 0.2
                    continue
                end

                V = 1.0im .* λ(R0[:,i]) .* (σX.*δR[2,i] .- σY.*δR[1,i])

                for i0=1:d, j0=1:d
                    if abs(V[i0,j0]) < tmin
                        continue
                    end
                    count = count+1
                    IS[count], JS[count], VS[count] = (i-1)*d + i0,  (j-1)*d + j0, V[i0,j0]
                end
            end
        end

        @views sparse(IS[1:count],JS[1:count],complex(VS[1:count]), d*N, d*N)
    end

    Hops(lat, hop; vectorized=true, cellrange=cellrange, format=format)
end
getrashba(lat::Lattice, λ::AbstractFloat; kwargs...) = getrashba(lat, x->λ; kwargs...)


addrashba!(hops, lat, rashba::Function; kwargs...) = addhops!(hops, getrashba(lat, rashba; kwargs...))
function addrashba!(hops, lat, rashba::Number)
    if !isapprox(rashba, 0; atol=sqrt(eps()))
        t(x) = rashba
        addrashba!(hops, lat, t)
    end
end

# One hopping law for scalar and matrix construction paths.
@inline function _graphene_hopping(Δ, δz; tmin=1e-5, tz=0.46, t0=1.0,
    ℓinter=0.125, ℓintra=0.08, z=3.0, a=1.0, Δmin=0.1, Δmax=5.0)
    (Δ == 0 || Δ < Δmin || Δ > Δmax) && return 0.0
    χ = δz^2 / Δ^2
    v = -t0 * (1 - χ) * exp(-abs(Δ - a) / ℓintra) -
        tz * χ * exp(-abs(Δ - z) / ℓinter)
    abs(v) < tmin ? zero(v) : v
end

function t_graphene(r1::AbstractVector, r2::AbstractVector; kwargs...)
    δx, δy, δz = r1[1] - r2[1], r1[2] - r2[2], r1[3] - r2[3]
    _graphene_hopping(sqrt(abs2(δx) + abs2(δy) + abs2(δz)), δz; kwargs...)
end

function t_graphene(R1::AbstractMatrix, R2::AbstractMatrix; kwargs...)
    N = size(R1, 2)
    out = N <= MAX_DENSE ? zeros(ComplexF64, N, N) : spzeros(ComplexF64, N, N)
    N > MAX_DENSE && sizehint!(out, MAX_DIAGS * N)
    addgraphenematrix!(out, R1, R2; kwargs...)
end

function addgraphenematrix!(out::AbstractMatrix, R1::AbstractMatrix, R2::AbstractMatrix; kwargs...)
    N = size(R1, 2)
    @assert N == size(out, 1) == size(out, 2) == size(R2, 2) "Incompatible dimensions."
    for j in 1:N, i in 1:N
        v = t_graphene(view(R1, :, i), view(R2, :, j); kwargs...)
        iszero(v) || (out[i, j] += v)
    end
    out
end
