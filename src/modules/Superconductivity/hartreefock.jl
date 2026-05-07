import LatticeQM.Meanfield
import LatticeQM.Utils
import SparseArrays: SparseMatrixCSC, findnz

mutable struct HartreeFockBDG{K, T2h, Th<:Hops{K,T2h}, T2v, Tv<:Hops{K,T2v}} <: Meanfield.MeanfieldGenerator{Th}
    const h::Th
    const v::Tv
    μ::Float64
    const V0::T2v
    const hMF::Th
    const ΔMF::Th
    ϵH::Float64
    ϵF::Float64
    ϵP::Float64
    ϵband::Float64
    ϵkin::Float64
    const fock::Bool
    const hartree::Bool

    function HartreeFockBDG(h::Th, v::Tv, μ=0.0; hartree=true, fock=true) where {K,T2h,Th<:Hops{K,T2h},T2v,Tv<:Hops{K,T2v}}
        V0 = sum(v[L] for L in keys(v))
        # See note in `Meanfield.HartreeFock`: type-preserving zero through
        # the underlying matrix, since `Base.zero(::Hops)` always densifies.
        hMF = Th(Dict{K,T2h}(L => zero(h[L]) for L in keys(h)))
        ΔMF = Th(Dict{K,T2h}(L => zero(h[L]) for L in keys(h)))
        # `ϵkin` is set by the SCF driver and shares the BdG quasi-particle
        # convention of the underlying band energy; for BdG the variational
        # identity has extra pairing-channel pieces, so interpret with care.
        new{K,T2h,Th,T2v,Tv}(h, v, μ, V0, hMF, ΔMF, 0.0, 0.0, 0.0, 0.0, 0.0, fock, hartree)
    end
end

Meanfield.sanitize!(ρ::BdGOperator{T}) where {T} = Meanfield.sanitize!(ρ.h)

# Anderson acceleration for BdG SCF: the fixedpoint! driver projects the
# post-Anderson iterate back onto the Hermitian subspace via
# `_hermitianize_if_supported!`. The base Meanfield module knows how to do
# this for plain `Hops` but not for `BdGOperator` (which lives downstream
# here), so without this method the projection silently no-ops on Nambu
# density matrices.
#
# Why the projection is needed (it's not paranoia): for Hermitian-symmetric
# iterates the Anderson coefficients γ from `ΔG \ g_k` are mathematically
# real, but Julia's complex `\` doesn't exploit that structure and returns γ
# with O(ε_machine) imaginary part. That seeds an O(ε) non-Hermitian
# component into x_{k+1}, the history then contains slightly non-Hermitian
# iterates, and the next γ picks up a larger imaginary part from
# `Re(ΔG^H ΔG)` no longer being purely real. Empirically the Hermiticity
# residue grows by ~10× per iteration on the BdG SCF — after ~17 steps it
# crosses `sqrt(eps()) ≈ 1.5e-8` and the next `BdGOperator(hMF, ΔMF)`
# constructor's hermiticity assertion fires. Projecting onto the Hermitian
# manifold each step holds the residue at machine precision indefinitely.
#
# A more principled alternative would be to enforce real γ at the LS-solve
# level (split into real/imag and solve a 2× real system); the projection
# achieves the same end state more cheaply.
Meanfield._hermitianize_if_supported!(x::BdGOperator) =
    (TightBinding.hermitianize!(x.h); x)

function Meanfield.HartreeFock(h::BdGOperator{Th}, v::Tv, μ=0.0; hartree=true, fock=true) where {K,T2h,Th<:Hops{K,T2h},T2v,Tv<:Hops{K,T2v}}
    h_el = Utils.copyelectronsector(h)
    HartreeFockBDG(h_el, v, μ; hartree=hartree, fock=fock)
end

function Meanfield.hMF(hf::HartreeFockBDG)
    hBdG = BdGOperator(hf.hMF, hf.ΔMF)
    if Meanfield.HARTREEFOCK_DEBUG[]
        @assert TightBinding.ishermitian(hBdG) "Mean-field Hamiltonian is not hermitian"
    end
    hBdG
end

function (hf::HartreeFockBDG)(ρ)
    if Meanfield.HARTREEFOCK_DEBUG[]
        @assert TightBinding.ishermitian(ρ) "Density matrix ρ is not hermitian"
    end
    ρ_el = getelectronview(ρ)
    ρ_Δ = getpairingview(ρ)

    meanfieldOperator!(hf, ρ_el, ρ_Δ)
    meanfieldScalar!(hf, ρ_el, ρ_Δ)
    hf
end

function initialize_ΔMF!(hf::HartreeFockBDG, ρΔ)
    T = eltype(hf.h[zerokey(hf.h)])
    dim = size(hf.h[zerokey(hf.h)])
    for k in keys(hf.ΔMF)
        hf.ΔMF[k] .= 0
    end
    for k in setdiff(union(keys(hf.h), keys(hf.v), keys(ρΔ)), keys(hf.ΔMF))
        hf.ΔMF[k] = zeros(T, dim)
    end
    nothing
end

function meanfieldOperator!(hf::HartreeFockBDG, ρ_el, ρΔ)
    Meanfield.initialize_hMF!(hf, ρ_el)
    initialize_ΔMF!(hf, ρΔ)

    if hf.fock
        Meanfield.meanfieldOperator_addfock!(hf, ρ_el)
        meanfieldOperator_addfock_pairing!(hf, ρΔ)
    end
    if hf.hartree
        Meanfield.meanfieldOperator_addhartree!(hf, ρ_el)
    end
    nothing
end

function meanfieldScalar!(hf::HartreeFockBDG, ρ_el, ρΔ)
    hf.ϵH = hf.hartree ? Meanfield.hartree_energy(hf, ρ_el) : 0.0
    hf.ϵF = hf.fock    ? Meanfield.fock_energy(hf, ρ_el)    : 0.0
    # Pairing channel: `meanfieldScalar_fock_pairing` returns `-½ Σ |Δ|² v`
    # (the double-counting form), so the physical pairing energy is its
    # negative.
    hf.ϵP = hf.fock    ? -meanfieldScalar_fock_pairing(hf, ρΔ) : 0.0
    nothing
end

function meanfieldOperator_addfock_pairing!(hf::HartreeFockBDG, ρΔ)
    for L in keys(hf.v)
        vL = hf.v[L]
        ΔL = ρΔ[L]
        ΔMF_L = hf.ΔMF[L]

        axes(vL) == axes(ΔL) || throw(DimensionMismatch("Pairing Fock block axis mismatch for key $(L): axes(v)=$(axes(vL)) vs axes(Δ)=$(axes(ΔL))"))
        axes(vL) == axes(ΔMF_L) || throw(DimensionMismatch("Pairing Fock block axis mismatch for key $(L): axes(v)=$(axes(vL)) vs axes(ΔMF)=$(axes(ΔMF_L))"))

        if vL isa SparseMatrixCSC
            rows, cols, vals = findnz(vL)
            @inbounds for idx in eachindex(vals)
                i = rows[idx]
                j = cols[idx]
                ΔMF_L[i, j] += vals[idx] * ΔL[i, j]
            end
        else
            @. ΔMF_L += vL * ΔL
        end
    end
    nothing
end

function meanfieldScalar_fock_pairing(hf::HartreeFockBDG, ρΔ)
    energy = -1/2 * sum(sum(ρΔ[L] .* conj.(ρΔ[L]) .* vL for (L, vL) in hf.v))
    @assert isapprox(imag(energy), 0; atol=Meanfield._imag_tol(ρΔ))
    real(energy)
end
