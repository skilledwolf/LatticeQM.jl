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
    ExchangeKernel(j, siteof, flavorof)

Inter-site exchange integrals for J-augmented Hartree-Fock. `j` is a
`Hops`-shaped kernel over **site** indices (dimension = number of sites per
cell, NOT the full orbital⊗flavor dimension), with `j[L][a,b]` the exchange
integral `J_ab(L) = ⟨a,0; b,L|V|b,L; a,0⟩ ≥ 0` between site `a` in the home
cell and site `b` in cell `L`. The self term `j[0][a,a]` must be zero (it is
part of the direct interaction). `siteof[I]`/`flavorof[I]` map each internal
orbital index `I` of the Hamiltonian basis to its site and flavor; every
(site, flavor) combination must occur exactly once, and flavor labels must be
consistent across sites (the interaction is flavor-independent, so exchange
pairs equal flavor labels on the two sites).

The corresponding interaction term is
`H_J = ½ Σ_L Σ_ab J_ab(L) Σ_αβ c†_{aα,0} c†_{bβ,L} c_{aβ,0} c_{bα,L}`,
i.e. `-Σ J (2 S_a·S_b + n_a n_b/2)` for spin-1/2 flavors — the standard
ferromagnetic direct exchange between orthogonal orbitals.
"""
struct ExchangeKernel{K,T2j,Tj<:Hops{K,T2j}}
    j::Tj
    sf2idx::Matrix{Int}   # (site, flavor) -> internal orbital index
    Jtot::T2j             # Σ_L j[L]
end

function ExchangeKernel(j::Hops, siteof::AbstractVector{<:Integer},
                        flavorof::AbstractVector{<:Integer})
    length(siteof) == length(flavorof) ||
        throw(ArgumentError("siteof and flavorof must have equal length"))
    ns, nf = maximum(siteof), maximum(flavorof)
    zk = zerokey(j)
    size(j[zk], 1) == ns ||
        throw(ArgumentError("exchange kernel dimension $(size(j[zk],1)) " *
                            "does not match number of sites $ns"))
    all(iszero, diag(j[zk])) ||
        throw(ArgumentError("onsite self-exchange j[0][a,a] must be zero " *
                            "(it is part of the direct interaction)"))
    sf2idx = zeros(Int, ns, nf)
    for I in eachindex(siteof)
        sf2idx[siteof[I], flavorof[I]] == 0 ||
            throw(ArgumentError("duplicate (site,flavor) = " *
                                "($(siteof[I]),$(flavorof[I]))"))
        sf2idx[siteof[I], flavorof[I]] = I
    end
    all(>(0), sf2idx) ||
        throw(ArgumentError("every (site,flavor) pair must appear once"))
    ExchangeKernel(j, sf2idx, sum(j[L] for L in keys(j)))
end

"""
    HartreeFock(h, v, μ=0.0; hartree=true, fock=true, exchange=nothing)

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

# Exchange augmentation (J-augmented HF)

Passing `exchange = ExchangeKernel(j, siteof, flavorof)` adds the inter-site
exchange interaction (see `ExchangeKernel`) at the same mean-field level.
Two fields are added to `hMF`, exact conjugates of the energy

  * `ϵJ = ½ Σ_L Σ_ab J_ab(L) ( |τ_L[a,b]|² − Tr_f[m_a m_b] )`

with `m_a` the onsite flavor matrix `ρ[0]` at site `a` and
`τ_L[a,b] = Σ_α ρ[L][(aα),(bα)]` the flavor-traced bond:

  * local flavor field  `Σ_loc(a) = −Σ_b (Σ_L J_ab(L)) m_b`  (matrix-valued
    Hartree; ferromagnetic — favors equal flavor polarization on all sites)
  * bond field  `Σ_bond[L][(aα),(bα)] = +J_ab(L) τ_L[a,b]`

The variational identity `Tr[Σ_J·ρ] = 2ϵJ` holds, all three total-energy
expressions above remain valid with `ϵH + ϵF → ϵH + ϵF + ϵJ`
(see `doublecounting`).
"""
mutable struct HartreeFock{K, T2h, Th<:Hops{K,T2h}, T2v, Tv<:Hops{K,T2v}, TJ<:Union{Nothing,ExchangeKernel}} <: MeanfieldGenerator{Th}
    const h::Th
    const v::Tv
    μ::Float64
    const V0::T2v
    const hMF::Th
    ϵH::Float64
    ϵF::Float64
    ϵJ::Float64
    ϵband::Float64
    ϵkin::Float64
    const fock::Bool
    const hartree::Bool
    const exchange::TJ

    function HartreeFock(h::Th, v::Tv, μ=0.0; hartree=true, fock=true, exchange=nothing) where {K,T2h,Th<:Hops{K,T2h},T2v,Tv<:Hops{K,T2v}}
        V0 = sum(v[L] for L in keys(v))
        if exchange !== nothing
            length(exchange.sf2idx) == size(h[zerokey(h)], 1) ||
                throw(ArgumentError("exchange site×flavor count " *
                    "$(length(exchange.sf2idx)) does not match Hamiltonian " *
                    "dimension $(size(h[zerokey(h)], 1))"))
        end
        # Pre-allocate hMF preserving h's matrix type — `Base.zero(::Hops)` is
        # overridden in Operators/densitymatrix.jl to always return dense
        # (because density-matrix partials must be dense), so we go through
        # the underlying matrix's `zero` instead to keep sparse Hamiltonians
        # sparse.
        hMF = Th(Dict{K,T2h}(L => zero(h[L]) for L in keys(h)))
        new{K,T2h,Th,T2v,Tv,typeof(exchange)}(h, v, μ, V0, hMF, 0.0, 0.0, 0.0,
                                              0.0, 0.0, fock, hartree, exchange)
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
    if hf.exchange !== nothing
        meanfieldOperator_addexchange!(hf, ρ)
    end

    nothing
end

function meanfieldScalar!(hf::HartreeFock, ρs)
    hf.ϵH = hf.hartree ? hartree_energy(hf, ρs) : 0.0
    hf.ϵF = hf.fock    ? fock_energy(hf, ρs)    : 0.0
    hf.ϵJ = hf.exchange !== nothing ? exchange_energy(hf, ρs) : 0.0
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
doublecounting(hf::HartreeFock) = hf.ϵH + hf.ϵF + hf.ϵJ


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

"""
    meanfieldOperator_addexchange!(hf, ρ)

Add the two exchange fields of the J-augmented functional (see the
`HartreeFock` docstring): the local flavor-matrix field
`Σ_loc(a) = −Σ_b Jtot[a,b] m_b` on the zero-key block and the flavor-traced
bond field `Σ_bond[L][(aα),(bα)] = +J[L][a,b] τ_L[a,b]`.
"""
function meanfieldOperator_addexchange!(hf, ρ)
    ek = hf.exchange
    ns, nf = size(ek.sf2idx)
    ix = ek.sf2idx
    zk = zerokey(ρ)
    ρ0 = ρ[zk]
    H0 = hf.hMF[zk]

    # local flavor-matrix field (ferromagnetic, matrix-valued Hartree)
    @inbounds for a in 1:ns, b in 1:ns
        Jab = ek.Jtot[a, b]
        iszero(Jab) && continue
        for α in 1:nf, β in 1:nf
            H0[ix[a, α], ix[a, β]] -= Jab * ρ0[ix[b, α], ix[b, β]]
        end
    end

    # flavor-traced bond field
    for L in keys(ek.j)
        haskey(ρ, L) || continue
        ρL = ρ[L]
        if !haskey(hf.hMF, L)
            hf.hMF[L] = zero(hf.h[zerokey(hf.h)])
        end
        hL = hf.hMF[L]
        jL = ek.j[L]
        @inbounds for a in 1:ns, b in 1:ns
            Jab = jL[a, b]
            iszero(Jab) && continue
            τ = zero(eltype(ρL))
            for α in 1:nf
                τ += ρL[ix[a, α], ix[b, α]]
            end
            for α in 1:nf
                hL[ix[a, α], ix[b, α]] += Jab * τ
            end
        end
    end
    nothing
end

"""
    exchange_energy(hf, ρ) → Real

Physical exchange energy of the J-augmented functional,
`ϵJ = ½ Σ_L Σ_ab J_ab(L) ( |τ_L[a,b]|² − Tr_f[m_a m_b] )`
(negative for polarized states; the local term wins for repulsive `J`).
"""
function exchange_energy(hf, ρs)
    ek = hf.exchange
    ns, nf = size(ek.sf2idx)
    ix = ek.sf2idx
    zk = zerokey(ρs)
    ρ0 = ρs[zk]

    energy = zero(ComplexF64)
    @inbounds for a in 1:ns, b in 1:ns
        Jab = ek.Jtot[a, b]
        iszero(Jab) && continue
        s = zero(ComplexF64)
        for α in 1:nf, β in 1:nf
            s += ρ0[ix[a, α], ix[a, β]] * ρ0[ix[b, β], ix[b, α]]
        end
        energy -= Jab * s / 2                       # −½ J Tr_f[m_a m_b]
    end
    for L in keys(ek.j)
        haskey(ρs, L) || continue
        ρL = ρs[L]
        jL = ek.j[L]
        @inbounds for a in 1:ns, b in 1:ns
            Jab = jL[a, b]
            iszero(Jab) && continue
            τ = zero(ComplexF64)
            for α in 1:nf
                τ += ρL[ix[a, α], ix[b, α]]
            end
            energy += Jab * abs2(τ) / 2             # +½ J |τ_L[a,b]|²
        end
    end
    @assert isapprox(imag(energy), 0; atol=_imag_tol(ρs))
    real(energy)
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
