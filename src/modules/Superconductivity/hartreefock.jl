import ..Meanfield
using ..Meanfield: hartreefock, hartreefock_pairing

import ..Spectrum

function Meanfield.hartreefock(H::BdGOperator{T}, v::AbstractHops) where T<:AbstractHops

    vMF, ΔMF, ϵMF = hartreefock_pairing(v)

    function H_op(ρ::BdGOperator{T}) where T<:AbstractHops

        ρ0 = Spectrum.getelectronsector(ρ)
        ρΔ = getpairingsector(ρ)

        MFOP = BdGOperator(vMF(ρ0), ΔMF(ρΔ))
        # MFOP = BdGOperator(vMF(ρ0))

        addhops(H,MFOP)
    end

    function E_scalar(ρ::BdGOperator{T}) where T<:AbstractHops
        ρ0 = Spectrum.getelectronsector(ρ)
        ρΔ = getpairingsector(ρ)

        ϵMF(ρ0, ρΔ)
    end

    H_op, E_scalar
end

# precompile(Meanfield.solvehartreefock, (BdGOperator{Hops{Matrix{ComplexF64}}}, BdGOperator{Hops{Matrix{ComplexF64}}}, BdGOperator{Hops{Matrix{ComplexF64}}}, Float64))