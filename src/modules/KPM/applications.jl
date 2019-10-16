function dos(H::T, n::Int, R::Int; kwargs...) where {T <: AbstractMatrix}
"""
    Compute a Chebyshev expansion of the density of states of
    a given Hamiltonian to order n by using the KPM with stochastic
    trace. The number of random vectors is R.

    [Eq. (107), Weisse et al, Rev. Mod. Phys. 78 275]
"""

    H1 = deepcopy(H)
    a,b = tounitrange!(H1)

    μ = get_expansion_stochastic(H1, n, R; kwargs...) ./ size(H1,2)

    μ, a, b
end

function ldos(H::T, site::Int, n::Int; kwargs...) where {T <: AbstractMatrix}
"""
    Compute a Chebyshev expansion of the local density of states
    of a given Hamiltonian at site 'site' to order n by using the
    KPM.

    [Eq. (113), Weisse et al, Rev. Mod. Phys. 78 275]
"""
    H1 = deepcopy(H)
    D = size(H1,2)
    a,b = tounitrange!(H1)

    α = zeros(ComplexF64, D)
    α[site] = 1.0

    get_expansion(α, H, n) ./ D,  a, b
end
