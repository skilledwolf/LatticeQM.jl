"""
    This method takes the Hamiltonian single-particle operator h and an
    interaction potential v and returns meanfield functionals
        ℋ, E  s.t.  h_mf = ℋ[ρ]  and  ϵ_scalar = E[ρ].

    These functionals can be used to search for a self-consistent solution
    using solve_selfconsistent(...).
"""

# function hartreefock(h::Function, v::AnyHops)
#     mf, E = hartreefock_k(v)
#     ℋ(ρ) = k -> (h(k) .+ mf(ρ)(k))

#     ℋ, E
# end

# function hartreefock_k(v::AnyHops)
#     vMF, ϵMF = hartreefock(v)
#     getbloch(vMF), ϵMF
# end

function hartreefock(h::AnyHops, v::AnyHops)
    vMF, ϵMF = hartreefock(v)

    hMF(ρ::AnyHops) = h + vMF(ρ)

    hMF, ϵMF
end


function hartreefock(v::AnyHops)
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ] that builds the mean field hamiltonian

        This may look harmless but requires a careful derivation.
    """

    V0 = sum(v[L] for L in keys(v))
    vmf = empty(v)

    function vMF(ρ::AnyHops)
        empty!(vmf)

        for L in keys(v)
            vmf[L] = -v[L] .* conj(ρ[L])#ρ[L] #conj(ρ[L]) #transpose(ρ[L]) # Fock contribution
        end

        addhops!(vmf, Hops([0,0] => spdiagm(0 => V0 * diag(ρ[[0,0]])))) # Hartree contribution

        vmf
    end

    function ϵMF(ρs::AnyHops)
        vρ = diag(ρs[[0,0]])

        energy = - 1/2 * (transpose(vρ) * V0 * vρ) # Hartree contribution
        energy +=  1/2 * sum(sum(ρs[L] .* conj.(ρs[L]) .* vL for (L,vL) in v)) # Fock contribution


        @assert isapprox(imag(energy),0; atol=sqrt(eps()))
        real(energy)
    end

    vMF, ϵMF
end


function hartreefock_pairing(v::AnyHops)
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ] that builds the mean field hamiltonian

        This may look harmless but requires a careful derivation.
    """

    V0 = sum(v[L] for L in keys(v))
    vmf = empty(v)
    Δmf = empty(v)

    function vMF(ρ::AnyHops)
        empty!(vmf)

        for L in keys(v)
            vmf[L] = -v[L] .* conj.(ρ[L])#ρ[L] #conj(ρ[L]) #transpose(ρ[L]) # Fock contribution
        end

        addhops!(vmf, Hops([0,0] => spdiagm(0 => V0 * diag(ρ[[0,0]])))) # Hartree contribution

        vmf
    end

    function ΔMF(ρ::AnyHops)
        empty!(Δmf)

        for L in keys(v)
            Δmf[L] = v[L] .* transpose(ρ[L])#ρ[L] #conj(ρ[L]) #transpose(ρ[L]) # Fock contribution
        end

        # addhops!(vmf, Hops([0,0] => spdiagm(0 => V0 * diag(ρ[[0,0]])))) # Hartree contribution

        Δmf
    end

    function ϵMF(ρs::AnyHops, ρΔs::AnyHops)
        vρ = diag(ρs[[0,0]])

        energy = - 1/2 * (transpose(vρ) * V0 * vρ) # Hartree contribution
        energy +=  1/2 * sum(sum(ρs[L] .* conj.(ρs[L]) .* vL for (L,vL) in v)) # Fock contribution
        energy -=  1/2 * sum(sum(ρΔs[L] .* conj.(ρΔs[L]) .* vL for (L,vL) in v)) # pairing contribution

        @assert isapprox(imag(energy),0; atol=sqrt(eps()))
        real(energy)
    end

    vMF, ΔMF, ϵMF
end