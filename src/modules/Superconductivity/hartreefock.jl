import ..Meanfield
# using ..Meanfield: hartreefock_pairing

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




function hartreefock_pairing(v::Hops)
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ] that builds the mean field hamiltonian

        This may look harmless but requires a careful derivation.
    """

    V0 = sum(v[L] for L in keys(v))
    # vmf = empty(v)
    # Δmf = empty(v)
    vmf = zero(v)
    Δmf = zero(v)

    function vMF(ρ::Hops)
        # empty!(vmf)
        fill!(vmf, 0.0)

        for L in keys(v)
            vmf[L] .+= -v[L] .* conj.(ρ[L]) # Fock contribution
        end

        vmf[zerokey(ρ)] += spdiagm(0 => V0 * diag(ρ[zerokey(ρ)])) # Hartree contribution
        # addhops!(vmf, Hops(zerokey(ρ) => spdiagm(0 => V0 * diag(ρ[zerokey(ρ)])))) # Hartree contribution

        vmf
    end

    function ΔMF(ρ::Hops)
        # empty!(Δmf)
        fill!(Δmf, 0.0)

        for L in keys(v)
            Δmf[L] .+= v[L] .* conj.(ρ[L]) # Fock contribution
        end

        Δmf
    end

    function ϵMF(ρs::Hops, ρΔs::Hops)
        vρ = diag(ρs[[0, 0]])

        energy = -1 / 2 * (transpose(vρ) * V0 * vρ) # Hartree contribution
        energy += 1 / 2 * sum(sum(ρs[L] .* conj.(ρs[L]) .* vL for (L, vL) in v)) # Fock contribution
        energy -= 1 / 2 * sum(sum(ρΔs[L] .* conj.(ρΔs[L]) .* vL for (L, vL) in v)) # pairing contribution

        @assert isapprox(imag(energy), 0; atol=sqrt(eps()))
        real(energy)
    end

    vMF, ΔMF, ϵMF
end