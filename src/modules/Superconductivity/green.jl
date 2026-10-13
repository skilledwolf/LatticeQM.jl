
import LatticeQM.Operators
import LatticeQM.TightBinding
import LatticeQM.Structure
import LinearAlgebra: tr
# using ..Operators: densitymatrix, densitymatrix!
# using ..TightBinding: zerokey
import SparseArrays

"""
    getdensitymatrix!(ρ::BdGOperator, H::BdGOperator, ks, μ=0; kwargs...)

Compute the full Nambu density matrix `ρ` from the BdG Hamiltonian `H` at
chemical potential `μ`. Internally shifts `H` by `−μ` while running the
sum, then restores `H` to its prior state (`try`/`finally`-protected).

Returns the physical canonical band energy `⟨H_phys⟩`, matching the normal
Hamiltonian density-matrix routine. The shifted Nambu sum is corrected by
its c-number and by `μ Nₑ`:

    ⟨H_phys⟩ = Tr(ρ_BdG H_BdG^shifted)/2 + (Tr h_e − N μ)/2 + μ Nₑ,

where `Nₑ = Tr(ρ_e[0])` is the physical electron number per cell.
"""
function Operators.getdensitymatrix!(ρ::BdGOperator{T}, H::BdGOperator{T2}, ks::AbstractMatrix, μ::Real=0; kwargs...) where {T,T2}
    L = size(ks, 2)
    kweights = fill(1 / L, L)
    Operators.getdensitymatrix!(ρ, H, ks, kweights, μ; kwargs...)
end

function Operators.getdensitymatrix!(ρ::BdGOperator{T}, H::BdGOperator{T2}, kgrid::Structure.Mesh, μ::Real=0; kwargs...) where {T,T2}
    Operators.getdensitymatrix!(ρ, H, kgrid.points, Structure.meshweights(kgrid), μ; kwargs...)
end

function Operators.getdensitymatrix!(ρ::BdGOperator{T}, H::BdGOperator{T2}, ks::AbstractMatrix, kweights::AbstractVector, μ::Real=0; kwargs...) where {T,T2}
    d = hopdim(H)

    Operators.addchemicalpotential!(H, -μ)
    out = try
        Operators.getdensitymatrix!(ρ.h, H.h, ks, kweights, 0.0; kwargs...)
    finally
        # Restore even on exception so caller's H is preserved.
        Operators.addchemicalpotential!(H, μ)
    end

    h_zk_electron = view(H.h[TightBinding.zerokey(H.h)], 1:d, 1:d)
    ρ_zk_electron = view(ρ.h[TightBinding.zerokey(ρ.h)], 1:d, 1:d)
    return out / 2 + (real(tr(h_zk_electron)) - d * μ) / 2 + μ * real(tr(ρ_zk_electron))
end
