import ..Meanfield
using ..Meanfield: hartreefock, hartreefock_pairing

function Meanfield.hartreefock(H::BdGOperator{T}, v::AnyHops) where T<:AnyHops

    vMF, ΔMF, ϵMF = hartreefock_pairing(v)

    function H_op(ρ::BdGOperator{T}) where T<:AnyHops

        ρ0 = getelectronsector(ρ)
        ρΔ = getpairingsector(ρ)

        MFOP = BdGOperator(vMF(ρ0), ΔMF(ρΔ))
        # MFOP = BdGOperator(vMF(ρ0))

        addhops(H,MFOP)
    end

    function E_scalar(ρ::BdGOperator{T}) where T<:AnyHops
        ρ0 = getelectronsector(ρ)
        ρΔ = getpairingsector(ρ)

        ϵMF(ρ0, ρΔ)
    end

    H_op, E_scalar
end