import LatticeQM.TightBinding
using LatticeQM.TightBinding: Hops, hopdim, addhops!, addhops, zerokey


"""
    BdGOperator(h::Hops)
    BdGOperator(h::Hops, Δ::Hops)

Bogoliubov–de Gennes operator wrapper that lifts a normal-state hopping object
`h` to Nambu space and optionally adds a pairing block `Δ`. The resulting
operator behaves like a Hamiltonian `H(k)` and can be passed to spectrum and
mean-field routines.

The Nambu Hamiltonian built is

    H_BdG(k) = [ h(k) − μ I    Δ(k)         ]
               [ Δ†(k)         −h(−k)^T + μ I ]

with the spinless / pre-spin-doubled convention `Ψ_i = (c_i, c_i^†)^T`. The
hole block is constructed as `−conj(h(R))`, which equals `−h(−R)^T` only
when `h` is Hermitian — the constructor asserts this. The chemical potential
is added later via `addchemicalpotential!(H::BdGOperator, μ)`, not at
construction time.

Notes on the dimensions:
- `hopdim(H::BdGOperator)` returns the **electron-sector** size N (the
  number of orbitals in the underlying `h`).
- `Spectrum.dim(H::BdGOperator)` returns the **Nambu** size 2N (what the
  eigensolver sees).

Related helpers:
- `addpairing!(H::BdGOperator, Δ)` — insert/update the pairing block.
- `addelectronsector!(H::BdGOperator, h)` — add to the electron block.
- `getelectronview` / `getpairingview` — views into electron/pairing sectors.
"""
mutable struct BdGOperator{T}
    h::T
end

(H::BdGOperator)(k) = H.h(k) # This will make the type callable, because T is callable


function BdGOperator(h0::T) where T<:Hops
    @assert TightBinding.ishermitian(h0) "BdGOperator(h): the normal-state hopping `h` must be Hermitian (h(R)† = h(−R)). Run TightBinding.hermitianize!(h) first if you have an asymmetric Hops."
    h1 = TightBinding.empty(h0) # create empty copy
    for R=keys(h0)
        # Lower-right block: -conj.(h(R)). Using h Hermitian, this equals
        # -h(-R)^T, which is the standard hole-block Fourier kernel that
        # produces -h^T(-k) in k-space.
        h1[R] = [h0[R]         zero(h0[R]);
                 zero(h0[R]) -conj.(h0[R])]
    end
    BdGOperator{T}(h1)
end

function BdGOperator(h0::T, Δ0::T2) where {T<:Hops,T2<:Hops}
    H = BdGOperator(h0)
    addpairing!(H, Δ0)
    H
end

"""
    addpairing!(H::BdGOperator, Δ::Hops)

Add or update the pairing block of a Nambu Hamiltonian. The user must
supply a `Δ` whose key set is closed under negation — every `R` in
`keys(Δ)` must have a partner `−R`. The block built at offset `R` is

    [ 0       Δ[ R] ]
    [ Δ[−R]† 0     ]

Note that `Δ[ R]` and `Δ[−R]` parametrise the *independent* upper-right
blocks at offsets `R` and `−R` (the lower-left at `R` is
`Δ[−R]†`). The resulting BdG is Hermitian as a Hops *regardless* of any
relationship between `Δ[ R]` and `Δ[−R]` — the only requirement is that
both keys are present so the construction can run. This deliberately
matches what the SCF produces: the mean-field pairing
`ΔMF[L] = v[L] .* ρΔ[L]` is generically *not* dict-Hermitian
(`ρΔ` is the upper-right block of `ρ_BdG`, which does not equal the
lower-left in general), but the BdG built from it is still Hermitian.
"""
function addpairing!(H::BdGOperator{T}, Δ::T2) where {T<:Hops, T2<:Hops}
    @assert hopdim(H) == hopdim(Δ) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    # Closure-under-negation check: addpairing! reads Δ[-R] when building
    # the lower-left block at R, so a missing -R key would crash with a
    # cryptic KeyError. Surface a clear message naming the offending offset.
    for R in keys(Δ)
        haskey(Δ, -R) || error("addpairing!: pairing dict has key $R but no partner $(- R). The construction reads Δ[-R] for the lower-left block at R — every R needs its −R counterpart present. Δ[ R] and Δ[−R] may take independent values; the BdG comes out Hermitian regardless.")
    end

    Δ1 = T2(Dict())
    for R=keys(Δ)
        Δ1[R] = [zero(Δ[R])   Δ[R];
                 Δ[-R]'       zero(Δ[R])]
    end

    addhops!(H,Δ1) # effectively mergewith!(+,)

    # Final guard: the construction is mathematically Hermitian when the
    # pre-existing h block is (asserted in BdGOperator(h)), so this only
    # fires on numerical drift or if a caller mutated H.h directly between
    # the BdGOperator(h) call and addpairing!.
    @assert TightBinding.ishermitian(H.h) "addpairing!: resulting BdG operator is not Hermitian. This indicates h drifted from Hermitian or H.h was mutated externally — the addpairing! construction itself preserves Hermiticity."

    H
end

function addelectronsector!(H::BdGOperator{T}, h::T2) where {T<:Hops,T2<:Hops}
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."
    @assert TightBinding.ishermitian(h) "addelectronsector!: the electron-sector hopping `h` must be Hermitian."

    H2 = T2(Dict())
    for R=keys(h)
        H2[R] = [h[R] zero(h[R]); zero(h[R]) -conj.(h[R])]
    end
    # H2 = Hops(Dict(R=>[h[R] zero(h[R]); zero(h[R]) -h[R]] for R=keys(h)))
    addhops!(H,H2)
end

function getelectronview(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)

    # hops = Hops()
    # for R = keys(H)
    #     hops[R] = view(H[R], 1:d, 1:d) #select electron part (no pairing, no holes)
    # end
    # hops
    Hops(Dict(R => view(H[R], 1:d, 1:d) for R = keys(H)))
end

function getpairingview(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)

    # Δ1 = Hops()
    # for R = keys(H)
    #     Δ1[R] = view(H[R], 1:d, d+1:2*d) #select electron part (no pairing, no holes)
    # end
    # Δ1
    Hops(Dict(R => view(H[R], 1:d, d+1:2*d) for R = keys(H)))
end


import LatticeQM.Utils

function Utils.getelectronsector(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)
    h = Dict(R=>view(H[R], 1:d, 1:d) for R=keys(H))
    Hops(h) # creates a view into the original object # TODO: check if this is works as intended
end
function Utils.copyelectronsector(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)
    Hops(Dict(R => H[R][1:d, 1:d] for R = keys(H))) # creates a copy of the original object # TODO: check if this is works as intended
end

function getpairingsector(H::BdGOperator{T}) where T<:Hops
    Δ1 = T(Dict())
    d = hopdim(H)

    for R=keys(H)
        Δ1[R] = H[R][1:d, d+1:2*d] #select electron part (no pairing, no holes)
    end

    Δ1
end


# Indexing interface
Base.values(H::BdGOperator) = Base.values(H.h)
Base.keys(H::BdGOperator) = Base.keys(H.h)
Base.getindex(H::BdGOperator,i) = Base.getindex(H.h,i)
Base.setindex!(H::BdGOperator,v,i) = Base.setindex!(H.h,v,i)
Base.firstindex(H::BdGOperator) = Base.firstindex(H.h)
Base.lastindex(H::BdGOperator) = Base.lastindex(H.h)


# Interface to tight binding module
TightBinding.hopdim(H::BdGOperator) = div(TightBinding.hopdim(H.h),2)

TightBinding.hermitianize!(H::BdGOperator) = TightBinding.hermitianize!(H.h)
TightBinding.ishermitian(H::BdGOperator) = TightBinding.ishermitian(H.h)

import ..Utils
import SparseArrays 
import ..TightBinding
Utils.dense(H::BdGOperator) = (h0 = Utils.dense(H.h); BdGOperator{typeof(h0)}(h0))
Utils.densecopy(H::BdGOperator) = (h0 = Utils.densecopy(H.h); BdGOperator{typeof(h0)}(h0))
SparseArrays.sparse(H::BdGOperator) = (h0 = SparseArrays.sparse(H.h); BdGOperator{typeof(h0)}(h0))
TightBinding.shareddense(H::BdGOperator) = (h0 = TightBinding.shareddense(H.h); BdGOperator{typeof(h0)}(h0))

Base.:-(h1::BdGOperator, h2::BdGOperator) = h1.h - h2.h
Base.:-(h1::BdGOperator, h2::Hops) = h1.h - h2
Base.:+(h1::BdGOperator, h2::BdGOperator) = h1.h + h2.h
Base.:+(h1::BdGOperator, h2::Hops) = h1.h + h2
TightBinding.addhops!(H1::BdGOperator, H2::BdGOperator) = TightBinding.addhops!(H1, H2.h)
TightBinding.addhops!(H1::BdGOperator, h::Hops) = (TightBinding.addhops!(H1.h, h); H1)
    

TightBinding.addhops(H1::BdGOperator, H2::BdGOperator) = TightBinding.addhops(H1, H2.h)
function TightBinding.addhops(H1::BdGOperator, h::Hops)
    Hnew = deepcopy(H1)
    TightBinding.addhops!(Hnew.h, h)
    Hnew
end

TightBinding.zerokey(H::BdGOperator) = TightBinding.zerokey(H.h)
