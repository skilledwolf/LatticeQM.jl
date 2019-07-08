function build_BlochH(hoppings::Dict{Vector{Int}, T1}; mode=:nospin) where {T2<:Complex, T1<:AbstractMatrix{T2}} #; mode=:auto
"""
    Returns the Bloch Hamiltonian with the BZ mapped onto the unit square, i.e., k ∈ [0,1]×[0,1].

    The input is a dictionary with keys being the reciprocal lattice vector in integer coordinates and
    the values are hopping matrices.
"""

    if mode==:spinhalf
        for (δL, t) in hoppings
            hoppings[δL] = kron(σ0, t) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
        end
    end

    BlochPhase(k,δL)::ComplexF64  = exp(1.0im * 2 * π * ComplexF64(dot(k,δL)))

    function hamiltonian(k::AbstractVector{Float64})#T3 where {T4<:Real, T3<:AbstractVector{T4}}

        res = sum(SparseMatrixCSC(t .* BlochPhase(k, δL) ) for (δL,t) in hoppings)
        res .= 0.5 .* (SparseMatrixCSC(res) .+ SparseMatrixCSC(res'))

        return res # not sure: shall I compute this on-the-fly (as here) or store the result?
    end

    hamiltonian
end
