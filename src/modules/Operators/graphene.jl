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
    hops!(hops, lat, t_graphene; kwargs...)
end

# default parameters taken from PHYSICAL REVIEW B 82, 035409 (2010), Table 1
function graphene_rhombohedral(lat; spin=false, a=1.0, d=3.0, γ0=-3.16, γ1=0.502, γ2=-0.0171, γ3=-0.377, γ4=-0.099, kwargs...)

    function t(r1, r2=0.0)
        δr=r1.-r2
        
        if abs(norm(δr[1:2])-a)<0.01 && abs(δr[3]) < 0.01 # same layer, NN
            return γ0
        elseif abs(norm(δr[1:2]))<0.01 && abs(abs(δr[3])-d) < 0.01  ## NN-layer, vertical
            return γ1
        elseif abs(norm(δr[1:2]))<0.01 && abs(abs(δr[3])-2*d) < 0.01 ## NNN-layer, vertical
            return γ2
        elseif abs(norm(δr[1:2])-a)<0.01 && abs(abs(δr[3])-d) < 0.01 && abs(δr[4]) < 0.01 ## NN-layer, non-vertical, AA
            return γ4
        elseif abs(norm(δr[1:2])-2*a)<0.01 && abs(abs(δr[3])-d) < 0.01 && abs(δr[4]) > 0.01 ## NN-layer, non-vertical, AB
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
`hops` on lattice `lat`. Positive values raise A and lower B (by convention).

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

################################################################################
################################################################################
################################################################################

"""
    This is the hopping amplitude between to twisted layers as given by overlap
    intergrals between p_z orbitals.
"""

norm2(r) = dot(r,r)

using ..Utils: heaviside

# graphene_interlayer(δz, Δ; z, ℓinter) =  δz^2 /Δ^2 * exp(-(Δ-z)/ℓinter)
# graphene_intralayer(δz, Δ; a, ℓintra, ℓz) = (1 - δz^2 /(Δ^2)) * exp(-(Δ-a)/ℓintra -δz^2 /ℓz^2)

@inline function graphene_hopping(r::AbstractVector{Float64}, r0=0; tz::Float64=0.46, t0::Float64=1.0, ℓinter::Float64=0.125, ℓintra::Float64=0.08, z::Float64=3.0, a::Float64=1.0,
           Δmin::Float64=0.1, Δmax::Float64=5.0)
    r = r .- r0
    @views Δ = sqrt.(sum(abs2,r[1:3]))
    result = 0.0
    if Δmax > Δ > Δmin
        @views δz = r[3]
        χ = δz^2 /(Δ^2)
        # result +=  (-t0) * (1-χ) * exp(-(Δ-a)/ℓintra) * exp(-δz^2 /ℓz^2) - tz * χ * exp(-(Δ-z)/ℓinter)
        result +=  (-t0) * (1-χ) * exp(-(Δ-a)/ℓintra) - tz * χ * exp(-(Δ-z)/ℓinter)
    end
    result
end

using ..TightBinding: MAX_DENSE, MAX_DIAGS

switchlin(x::AbstractFloat) =  ifelse(x < 0, zero(x), x)

function t_graphene(R1::Matrix{Float64}, R2::Matrix{Float64}; kwargs...)
    N = size(R1,2)

    if N < MAX_DENSE + 1
        out = zeros(ComplexF64, N,N)
    else
        out = spzeros(ComplexF64, N,N)
        sizehint!(out, round(Int, MAX_DIAGS * N)) # preallocate memory: important for huge sparse matrices
    end

    addgraphenematrix!(out, R1, R2; kwargs...)
end
precompile(t_graphene, (Matrix{Float64}, Matrix{Float64}))

function addgraphenematrix!(out::AbstractMatrix, R1::Matrix{Float64}, R2::Matrix{Float64}; tmin::Float64=1e-5, tz::Float64=0.46, t0::Float64=1.0,
    ℓinter::Float64=0.125, ℓintra::Float64=0.08, z::Float64=3.0, a::Float64=1.0,
    Δmin::Float64=0.1, Δmax::Float64=5.0)

    N = size(R1, 2)
    @assert N == size(out,1) == size(out,2) == size(R2,2) "Incompatible dimensions."

    δR = similar(R1)
    @fastmath @inbounds for j = 1:N
        @views δR .= R1 .- R2[:, j]

        for i = 1:N
            @views Δ = sqrt(sum(abs2, δR[1:3, i])) # absolute value of distance vector

            if Δmax < Δ || Δ < Δmin # cutoff lengths
                continue
            end

            δz = δR[3, i]
            χ = δz^2 / (Δ^2)

            v = -t0 * (1 - χ) * exp(-abs(Δ - a) / ℓintra) - tz * χ * exp(-abs(Δ - z) / ℓinter)

            if abs(v) < tmin # small energy cutoff for hop amplitudes
                continue
            end

            out[i, j] += v
        end
    end

    out
end
precompile(addgraphene!, (Matrix{Float64}, Matrix{Float64}, Matrix{Float64}))
precompile(addgraphene!, (Matrix{ComplexF64}, Matrix{Float64}, Matrix{Float64}))
precompile(addgraphene!, (SparseMatrixCSC{ComplexF64,Int64}, Matrix{Float64}, Matrix{Float64}))
precompile(addgraphene!, (SparseMatrixCSC{Float64,Int64}, Matrix{Float64}, Matrix{Float64}))
