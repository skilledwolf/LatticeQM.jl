
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

Returns the **physical** grand-canonical kinetic energy
⟨H_phys − μ N̂⟩, recovered from the BdG sum via

    ⟨H_phys − μ N̂⟩ = (1/2) Tr(ρ_BdG · H_BdG^shifted) + (Tr h_e − N μ) / 2,

where the constant term is the Nambu c-number from anticommuting
`c c† = 1 − c†c` in the hole sector. Without it the reported value
would equal twice the kinetic plus a μ- and N-dependent offset; this
matters whenever `ϵ_GS = ϵ_kin + ϵ_MF` is consumed (e.g. comparing total
energies between SCF runs at different parameters).
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

    # Recover ⟨H_phys − μ N̂⟩ from the BdG sum. After the restore, the
    # zero-key electron block holds the user's h_e (no μ shift), so we
    # subtract `d * μ` explicitly to account for what was used inside.
    h_zk_electron = view(H.h[TightBinding.zerokey(H.h)], 1:d, 1:d)
    return out / 2 + (real(tr(h_zk_electron)) - d * μ) / 2
end
