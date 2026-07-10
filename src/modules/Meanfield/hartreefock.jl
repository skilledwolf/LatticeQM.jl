"""
    This method takes the Hamiltonian single-particle operator h and an
    interaction potential v and returns meanfield functionals
        ℋ, E  s.t.  h_mf = ℋ[ρ]  and  ϵ_scalar = E[ρ].

    These functionals can be used to search for a self-consistent solution
    using solve_selfconsistent(...).
"""

import LatticeQM.TightBinding: zerokey
import SparseArrays: SparseMatrixCSC, findnz

"""
    MeanfieldGenerator

Abstract supertype for mean‑field functionals that map a density matrix `ρ` to
an effective single‑particle Hamiltonian and scalar energy contributions. Concrete
implementations include `HartreeFock` and `HartreeFockBDG`.
"""
abstract type MeanfieldGenerator{T} end

"""
    HARTREEFOCK_DEBUG[]

Set to `true` to re-enable per-iteration Hermiticity checks on the mean-field
operator (off by default for performance: the check on a `Hops` runs O(N²)
work and `hMF(::HartreeFock)` is called multiple times per SCF step).
"""
const HARTREEFOCK_DEBUG = Ref(false)

"""
    HartreeFock(h, v, μ=0.0; hartree=true, fock=true)

Mean-field functional for density (Hartree) and exchange (Fock) channels built
from a base Hamiltonian `h` and interaction kernels `v`. Calling the struct on
`ρ` updates the effective mean-field operator `hMF` and the scalar energy
fields `ϵH`, `ϵF` (see *Energy decomposition* below). The SCF driver fills in
`ϵband` and `ϵkin` once a density-matrix iteration has run.

`h` and `v` may use different matrix backends (e.g. dense `h` with sparse
Hubbard `v`); they only need to share the lattice key type and dimensions.

# Convention

`v` and `h` (and the resulting density matrix `ρ`) must use the **same
orbital basis**. In particular, if your model has spin, the basis must
already include it (e.g. via `TightBinding.addspin`); the Hartree term reads
diagonal occupations from `ρ` and assumes those are the densities `v` couples
to. For a spinful Hubbard `U n_↑ n_↓` model, build `v` so its matrix elements
between spin-↑ and spin-↓ orbitals encode `U`, not the same-spin diagonal.

# Energy decomposition (physical-sign convention)

After `hf(ρ)`:

  * `ϵH = +½ nᵀ V₀ n = ½ Tr[V_H ρ]` — Hartree energy (positive for
    repulsive `V₀`)
  * `ϵF = -½ Σ_L Σ_{ij} v_L[i,j] |ρ_L[i,j]|² = ½ Tr[V_F ρ]` — Fock /
    exchange energy (negative for repulsive `v`)

After one SCF iteration the driver sets:

  * `ϵband = Tr[hMF · ρ]` — sum of occupied band energies
  * `ϵkin = Tr[h · ρ] = ϵband - 2(ϵH + ϵF)` — kinetic / bare-Hamiltonian
    expectation value (derived from the variational identity
    `E_band = ⟨T⟩ + 2(E_H + E_F)`)

The total HF energy can be read off in three equivalent ways — useful as
self-consistency checks at convergence:

  * `E_HF = ϵkin + ϵH + ϵF`            (physical, additive)
  * `E_HF = ϵband - ϵH - ϵF`           (band − double-counting)
  * `E_HF = ½ (ϵkin + ϵband)`          (`½ Tr[ρ (h₀ + hMF)]`)
"""
mutable struct HartreeFock{K, T2h, Th<:Hops{K,T2h}, T2v, Tv<:Hops{K,T2v}} <: MeanfieldGenerator{Th}
    const h::Th
    const v::Tv
    μ::Float64
    const V0::T2v
    const hMF::Th
    ϵH::Float64
    ϵF::Float64
    ϵband::Float64
    ϵkin::Float64
    const fock::Bool
    const hartree::Bool

    function HartreeFock(h::Th, v::Tv, μ=0.0; hartree=true, fock=true) where {K,T2h,Th<:Hops{K,T2h},T2v,Tv<:Hops{K,T2v}}
        V0 = sum(v[L] for L in keys(v))
        # Pre-allocate hMF preserving h's matrix type — `Base.zero(::Hops)` is
        # overridden in Operators/densitymatrix.jl to always return dense
        # (because density-matrix partials must be dense), so we go through
        # the underlying matrix's `zero` instead to keep sparse Hamiltonians
        # sparse.
        hMF = Th(Dict{K,T2h}(L => zero(h[L]) for L in keys(h)))
        new{K,T2h,Th,T2v,Tv}(h, v, μ, V0, hMF, 0.0, 0.0, 0.0, 0.0, fock, hartree)
    end
end

function hMF(hf::HartreeFock)
    if HARTREEFOCK_DEBUG[]
        @assert TightBinding.ishermitian(hf.hMF) "Mean-field Hamiltonian is not hermitian"
    end
    hf.hMF
end

function (hf::MeanfieldGenerator)(ρ)
    meanfieldOperator!(hf, ρ)
    meanfieldScalar!(hf, ρ)
    hf
end


function initialize_hMF!(hf, ρ)
    for k in keys(hf.hMF)
        hf.hMF[k] .= 0
    end
    # The Fock loop runs over keys(v); when the interaction reaches further
    # than the hopping (e.g. NN Hubbard V on an onsite-hopping model), hMF —
    # allocated with keys(h) — must be extended or the update KeyErrors.
    # (The BdG ΔMF channel already did this; the electron channel was missed.)
    for k in keys(hf.v)
        haskey(hf.hMF, k) && continue
        hf.hMF[k] = zero(hf.h[TightBinding.zerokey(hf.h)])
    end
    TightBinding.addhops!(hf.hMF, hf.h)
    nothing
end

function meanfieldOperator!(hf::HartreeFock, ρ)
    initialize_hMF!(hf, ρ)

    if hf.fock
        meanfieldOperator_addfock!(hf, ρ)
    end
    if hf.hartree
        meanfieldOperator_addhartree!(hf, ρ)
    end

    nothing
end

function meanfieldScalar!(hf::HartreeFock, ρs)
    hf.ϵH = hf.hartree ? hartree_energy(hf, ρs) : 0.0
    hf.ϵF = hf.fock    ? fock_energy(hf, ρs)    : 0.0
    nothing
end

"""
    doublecounting(hf) → Float64

Total interaction double-counting correction of the mean-field band energy:
`E_GS = ϵband − doublecounting(hf)` and `⟨T⟩ = ϵband − 2·doublecounting(hf)`.
For plain Hartree-Fock this is `ϵH + ϵF`; BdG generators add the pairing
channel `ϵP` (the quasiparticle band energy contains `⟨½(c†Δc† + h.c.)⟩ =
2ϵP` while the physical pairing interaction energy is `ϵP`).
"""
doublecounting(hf::MeanfieldGenerator) = hf.ϵH + hf.ϵF


####################################################################
# Low-level functions to construct Hartree and Fock contributions
####################################################################

function meanfieldOperator_addfock!(hf, ρ)
    for L in keys(hf.v)
        # ρ blocks the density matrix doesn't carry are identically zero, so
        # their exchange contribution vanishes — skip instead of KeyError.
        haskey(ρ, L) || continue
        vL = hf.v[L]
        ρL = ρ[L]
        hL = hf.hMF[L]

        axes(vL) == axes(ρL) || throw(DimensionMismatch("Fock block axis mismatch for key $(L): axes(v)=$(axes(vL)) vs axes(ρ)=$(axes(ρL))"))
        axes(vL) == axes(hL) || throw(DimensionMismatch("Fock block axis mismatch for key $(L): axes(v)=$(axes(vL)) vs axes(hMF)=$(axes(hL))"))

        # Exchange field conjugate to E_F = -½ Σ v|ρ|²: Σ_F[L] = -v_L .* ρ_L
        # (NOT conj(ρ_L) — both are Hermitian as Hops, but only the
        # un-conjugated form satisfies Tr[Σ_F·ρ] = 2ϵF; the conjugated one
        # iterates the time-reversed exchange field, so complex-ρ states
        # (spiral/XY order, loop currents) are not fixed points.)
        if vL isa SparseMatrixCSC
            rows, cols, vals = findnz(vL)
            @inbounds for idx in eachindex(vals)
                i = rows[idx]
                j = cols[idx]
                hL[i, j] += -vals[idx] * ρL[i, j]
            end
        else
            @. hL -= vL * ρL
        end
    end
    nothing
end

function meanfieldOperator_addhartree!(hf, ρ)
    # Add Hartree shift to the diagonal of the zero-key block. Direct
    # element-wise loop avoids constructing an intermediate sparse `spdiagm`
    # and the sparse-to-dense add when hMF is dense.
    zk = zerokey(ρ)
    H0 = hf.hMF[zk]
    d = hf.V0 * diag(ρ[zk])
    @inbounds for i in eachindex(d)
        H0[i, i] += d[i]
    end
    nothing
end

# Tolerance for "imag(energy) ≈ 0" assertions: the imaginary part is the
# residue of finite-precision arithmetic over O(N) accumulations, so we
# scale with the number of orbitals in the zero-key block.
_imag_tol(ρs) = sqrt(eps()) * max(1, size(ρs[zerokey(ρs)], 1))

"""
    hartree_energy(hf, ρ) → Real

Physical Hartree energy `+½ nᵀ V₀ n` (positive for repulsive `V₀`).
"""
function hartree_energy(hf, ρs)
    vρ = diag(ρs[zerokey(ρs)])
    energy = 1/2 * (transpose(vρ) * hf.V0 * vρ)
    @assert isapprox(imag(energy), 0; atol=_imag_tol(ρs))
    real(energy)
end

"""
    fock_energy(hf, ρ) → Real

Physical Fock (exchange) energy `-½ Σ_L Σ_{ij} v_L[i,j] |ρ_L[i,j]|²`
(negative for repulsive `v`).
"""
function fock_energy(hf, ρs)
    energy = -1/2 * sum(haskey(ρs, L) ? sum(ρs[L] .* conj.(ρs[L]) .* vL) : 0.0
                        for (L, vL) in hf.v)
    @assert isapprox(imag(energy), 0; atol=_imag_tol(ρs))
    real(energy)
end

# Legacy double-counting-correction helpers (negatives of the physical
# energies above). Kept because `Superconductivity.HartreeFockBDG` reuses
# them in its mean-field scalar; new code should call `hartree_energy` /
# `fock_energy` instead.
meanfieldScalar_hartree(hf, ρs) = -hartree_energy(hf, ρs)
meanfieldScalar_fock(hf, ρs)    = -fock_energy(hf, ρs)
