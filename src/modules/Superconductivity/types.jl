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

function BdGOperator(h0::T, Δ0::T) where {T<:Hops}
    H = BdGOperator(h0)
    addpairing!(H, Δ0)
    H
end

"""
    addpairing!(H::BdGOperator, Δ::Hops)

Add or update the pairing block of a Nambu Hamiltonian. The user is
responsible for providing a `Δ` that yields a Hermitian `H_BdG` — i.e. the
dictionary must contain `−R` for every `R` it contains, and
`Δ[−R] == Δ[R]^†` must hold. Both conditions are validated; mismatches raise
an error rather than being silently symmetrised, since the symmetrisation
choice depends on the physical model (s-wave / d-wave / spin-singlet vs
triplet decomposition).
"""
function addpairing!(H::BdGOperator{T}, Δ::T2) where {T<:Hops, T2<:Hops}
    @assert hopdim(H) == hopdim(Δ) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    # Validate the pairing dictionary up front, with messages that name the
    # offset(s) the user has to fix.
    for R in keys(Δ)
        haskey(Δ, -R) || error("addpairing!: pairing dict has key $R but no partner $(- R). For a Hermitian BdG, every R needs its −R counterpart with Δ[−R] = Δ[R]†. Add the missing offset (or use .−R = Δ[R]') before calling addpairing!.")
    end
    for R in keys(Δ)
        # Δ[-R] == Δ[R]† is required for Hermiticity of H_BdG. Check via a
        # max-abs threshold — strict isapprox would be too noisy at machine
        # epsilon for matrices built by user code.
        maxdev = maximum(abs, Δ[-R] .- Δ[R]')
        maxdev > sqrt(eps()) && error("addpairing!: Δ[−R] is not Δ[R]† for R=$R (max abs deviation $(maxdev)). The user must make the pairing dict Hermitian-compatible before calling addpairing!.")
    end

    Δ1 = T2(Dict())
    for R=keys(Δ)
        Δ1[R] = [zero(Δ[R])   Δ[R];
                 Δ[-R]'       zero(Δ[R])]
    end

    addhops!(H,Δ1) # effectively mergewith!(+,)

    # Final guard: assert that the resulting Nambu operator really is
    # Hermitian, in case the user-supplied h or Δ had numerical drift not
    # caught by ishermitian.
    @assert TightBinding.ishermitian(H.h) "addpairing!: resulting BdG operator is not Hermitian. This usually indicates Δ[−R] ≠ Δ[R]† at some offset, or h that drifted from Hermitian."

    H
end

function addelectronsector!(H::BdGOperator{T}, h::T2) where {T<:Hops,T2<:Hops}
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    H2 = T2(Dict())
    for R=keys(h)
        H2[R] = [h[R] zero(h[R]); zero(h[R]) -h[R]]
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

function TightBinding.gethopsview(H::BdGOperator)
    h_view = TightBinding.gethopsview(H.h)
    BdGOperator{typeof(h_view)}(h_view)
end

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
