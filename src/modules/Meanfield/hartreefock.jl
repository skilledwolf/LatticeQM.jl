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
    HartreeFock(h, v, μ=0.0; hartree=true, fock=true)

Mean‑field functional for density (Hartree) and exchange (Fock) channels built
from a base Hamiltonian `h` and interaction kernels `v`. Calling the struct on
`ρ` updates the effective mean‑field operator `hMF` and scalar energy `ϵMF`.
"""
mutable struct HartreeFock{K, T2, T<:Hops{K,T2}} <: MeanfieldGenerator{T} 
    h::T
    v::T
    μ::Float64
    V0::T2 
    hMF::T
    ϵMF::Float64
    fock::Bool
    hartree::Bool

    function HartreeFock(h::T, v::T, μ=0.0; hartree=true, fock=true) where {K,T2,T<:Hops{K,T2}}
        # Extract K and T2 from T if needed inside the constructor
        # Example: Suppose AbstractHops is defined as AbstractHops{K, T2}
        # K = Base.parameter_upper_bound(T, 1)
        # T2 = Base.parameter_upper_bound(T, 2)
        V0 = sum(v[L] for L in keys(v))
        hMF = zero(h)  # Pre-allocated for performance
        ϵMF = 0.0
        new{K,T2,T}(h, v, μ, V0, hMF, ϵMF, fock, hartree)
    end
end

function hMF(hf::HartreeFock)
    h0 = hf.hMF
    @assert TightBinding.ishermitian(h0) "Mean-field Hamiltonian is not hermitian"
    h0
    # hf.hMF
end

function (hf::MeanfieldGenerator)(ρ) #(ρ::T) where {K,T2,T<:Hops{K,T2}}
    meanfieldOperator!(hf, ρ)
    meanfieldScalar!(hf, ρ)
    hf
end


function initialize_hMF!(hf, ρ)

    for k in keys(hf.hMF)
        hf.hMF[k] .= 0 # reset hMF to zero
    end

    TightBinding.addhops!(hf.hMF, hf.h)

    # T = eltype(hf.h[zerokey(hf.h)])
    # dim = size(hf.h[zerokey(hf.h)])
    # for k in union(keys(hf.h), keys(hf.v), keys(ρ)) # initialize hMF with h
    #     if haskey(hf.hMF, k)
    #         if haskey(hf.h, k)
    #             hf.hMF[k] .= hf.h[k]
    #         else
    #             hf.hMF[k] .= 0
    #         end
    #     else
    #         if haskey(hf.h, k)
    #             hf.hMF[k] = hf.h[k]
    #         else
    #             # error("Implementation error: We should better not reach this part...")
    #             hf.hMF[k] = zeros(T, dim) # NOTE: for sparse operators, this will fail, but should not be reached usually?
    #         end
    #     end
    # end
    nothing 
end

function meanfieldOperator!(hf::HartreeFock, ρ)
    initialize_hMF!(hf, ρ) # Call the new initialization function

    if hf.fock
        meanfieldOperator_addfock!(hf, ρ)
    end
    if hf.hartree
        meanfieldOperator_addhartree!(hf, ρ)
    end

    # @todo: Potentially we should add a trimming stage here
    # for sparse matrices (i.e., drop entries numerically close to 0)
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
            hL .+= -vL .* conj.(ρL)
        end
    end
    nothing
end

function meanfieldOperator_addfock_pairing!(hf, ρΔ)
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
                ΔMF_L[i, j] += vals[idx] * conj(ΔL[i, j])
            end
        else
            ΔMF_L .+= vL .* conj.(ΔL)
        end
    end
    nothing
end

function meanfieldOperator_addhartree!(hf, ρ)
    hf.hMF[zerokey(ρ)] .+= spdiagm(0 => hf.V0 * diag(ρ[zerokey(ρ)])) # Hartree contribution
    nothing
end

function meanfieldScalar_hartree(hf, ρs)
    vρ = diag(ρs[zerokey(ρs)])
    energy = -1/2 * (transpose(vρ) * hf.V0 * vρ) # Hartree contribution
    @assert isapprox(imag(energy), 0; atol=sqrt(eps()))
    real(energy)
end

function meanfieldScalar_fock(hf, ρs)
    energy = 1/2 * sum(sum(ρs[L] .* conj.(ρs[L]) .* vL for (L, vL) in hf.v)) # Fock contribution
    @assert isapprox(imag(energy), 0; atol=sqrt(eps()))
    real(energy)
end

function meanfieldScalar_fock_pairing(hf, ρΔ)
    energy = -1/2 * sum(sum(ρΔ[L] .* conj.(ρΔ[L]) .* vL for (L, vL) in hf.v)) # Fock contribution
    @assert isapprox(imag(energy), 0; atol=sqrt(eps()))
    real(energy)
end
