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
`ρ` updates the effective mean-field operator `hMF` and scalar energy `ϵMF`.

`h` and `v` may use different matrix backends (e.g. dense `h` with sparse
Hubbard `v`); they only need to share the lattice key type and dimensions.

# Convention

`v` and `h` (and the resulting density matrix `ρ`) must use the **same
orbital basis**. In particular, if your model has spin, the basis must
already include it (e.g. via `TightBinding.addspin`); the Hartree term reads
diagonal occupations from `ρ` and assumes those are the densities `v` couples
to. For a spinful Hubbard `U n_↑ n_↓` model, build `v` so its matrix elements
between spin-↑ and spin-↓ orbitals encode `U`, not the same-spin diagonal.
"""
mutable struct HartreeFock{K, T2h, Th<:Hops{K,T2h}, T2v, Tv<:Hops{K,T2v}} <: MeanfieldGenerator{Th}
    const h::Th
    const v::Tv
    μ::Float64
    const V0::T2v
    const hMF::Th
    ϵMF::Float64
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
        ϵMF = 0.0
        new{K,T2h,Th,T2v,Tv}(h, v, μ, V0, hMF, ϵMF, fock, hartree)
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
    hf.ϵMF = 0.0

    if hf.fock
        hf.ϵMF += meanfieldScalar_fock(hf, ρs)
    end
    if hf.hartree
        hf.ϵMF += meanfieldScalar_hartree(hf, ρs)
    end
    nothing
end


####################################################################
# Low-level functions to construct Hartree and Fock contributions
####################################################################

function meanfieldOperator_addfock!(hf, ρ)
    for L in keys(hf.v)
        vL = hf.v[L]
        ρL = ρ[L]
        hL = hf.hMF[L]

        axes(vL) == axes(ρL) || throw(DimensionMismatch("Fock block axis mismatch for key $(L): axes(v)=$(axes(vL)) vs axes(ρ)=$(axes(ρL))"))
        axes(vL) == axes(hL) || throw(DimensionMismatch("Fock block axis mismatch for key $(L): axes(v)=$(axes(vL)) vs axes(hMF)=$(axes(hL))"))

        if vL isa SparseMatrixCSC
            rows, cols, vals = findnz(vL)
            @inbounds for idx in eachindex(vals)
                i = rows[idx]
                j = cols[idx]
                hL[i, j] += -vals[idx] * conj(ρL[i, j])
            end
        else
            @. hL -= vL * conj(ρL)
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

function meanfieldScalar_hartree(hf, ρs)
    vρ = diag(ρs[zerokey(ρs)])
    energy = -1/2 * (transpose(vρ) * hf.V0 * vρ) # Hartree contribution
    @assert isapprox(imag(energy), 0; atol=_imag_tol(ρs))
    real(energy)
end

function meanfieldScalar_fock(hf, ρs)
    energy = 1/2 * sum(sum(ρs[L] .* conj.(ρs[L]) .* vL for (L, vL) in hf.v)) # Fock contribution
    @assert isapprox(imag(energy), 0; atol=_imag_tol(ρs))
    real(energy)
end
