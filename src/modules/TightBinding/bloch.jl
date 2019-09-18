BlochPhase(k,δL)::ComplexF64  = exp(1.0im * 2 * π * ComplexF64(dot(k,δL)))


function get_bloch(hoppings::LatticeHops; mode=:nospin, symmetric=true) #; mode=:auto
"""
    Returns the Bloch Hamiltonian with the BZ mapped onto the unit square, i.e., k ∈ [0,1]×[0,1].

    The input is a dictionary with keys being the reciprocal lattice vector in integer coordinates and
    the values are hopping matrices.
"""

    # hoppings = extend_space(hoppings, mode)

    function hamiltonian(k::AbstractVector{Float64})#T3 where {T4<:Real, T3<:AbstractVector{T4}}

        res = sum(t .* BlochPhase(k, δL) for (δL,t) in hoppings)
        if symmetric
            res .= 0.5 .* (res .+ res')
        end

        return res # not sure: shall I compute this on-the-fly (as here) or store the result?
    end

    hamiltonian
end


function build_BlochH(args...; kwargs...)
    @warn("Deprecation warning: build_BlochH() was renamed to get_bloch() and is marked for removal.")
    get_bloch(args...; kwargs...)
end
