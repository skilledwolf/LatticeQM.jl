"""
    This method takes the Hamiltonian single-particle operator h and an
    interaction potential v and returns meanfield functionals
        ℋ, E  s.t.  h_mf = ℋ[ρ]  and  ϵ_scalar = E[ρ].

    These functionals can be used to search for a self-consistent solution
    using solve_selfconsistent(...).
"""
function hartreefock(h::Function, v::AnyHops)
    mf, E = hartreefock_k(v)
    ℋ(ρ) = k -> (h(k) .+ mf(ρ)(k))

    ℋ, E
end

hartreefock(h::AnyHops, v::AnyHops) = ρ::AnyHops -> addhops(h, hartreefock(v)(ρ))

function hartreefock_k(v::AnyHops)
    vMF, ϵMF = hartreefock(v)
    getbloch(vMF), ϵMF
end

function hartreefock(v::AnyHops)
    """
        Expects the real space potential {V(L) | L unit cell vector}.
        It returns a functional 𝒱[ρ] that builds the mean field hamiltonian

        This may look harmless but requires a careful derivation.
    """

    V0 = sum(v[L] for L in keys(v))

    function vMF(ρ::AnyHops)
        v = Hops([0,0] => spdiagm(0 => V0 * diag(ρ[[0,0]]))) # Hartree contribution

        for (L,ρL) in ρ
            v[L] .*= -v[L] .* ρ[L] # Fock contribution
        end
        v
    end

    function ϵMF(ρs::AnyHops)
        vρ = diag(ρs[[0,0]])

        hartree = - 1/2 * (transpose(vρ) * V0 * vρ) # Hartree contribution
        fock =  1/2 * sum(sum(ρL .* conj.(ρL) .* v[L] for (L,ρL) in ρs)) # Fock contribution

        @assert imag(hartree) ≈ 0 && imag(fock) ≈ 0
        real(hartree + fock)
    end

    vMF, ϵMF
end


###################################################################################################
# Backwards compatibility
###################################################################################################
@legacyalias hartreefock get_mf_functional
@legacyalias hartreefock_k get_mf_operator
