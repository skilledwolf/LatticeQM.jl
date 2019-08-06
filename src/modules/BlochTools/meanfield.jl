
using NLsolve

#AbstractArray{Float64,N}
function solve_selfconsistent(ℋ::Function,
    ρ::AbstractArray{T,N}, ks::AbstractMatrix{Float64}, filling::Float64;
    iterations=500, ftol=1e-7, xtol=1e-7, method=:anderson, m=5, format=:dense, kwargs...) where {N, T<:Number}

    type = issparse(ℋ(zero(ρ))(ks[:,1])) ? :sparse : :dense # Decide if the Hamiltonian is sparse

    function f!(ρ1, ρ0)
        # Update meanfield Hamiltonian and chemical potential
        h = ℋ(ρ0)

        Σ = spectrum(h; format=format)
        μ = chemical_potential(h, ks, filling; type=type)

        # Obtain the meanfield of the updated Hamiltonian
        density!(ρ1, Σ, ks, μ)#; format=format)

        # Update meanfield for next iteration step
        ρ1 .-= ρ0
        nothing
    end

    df = OnceDifferentiable(f!, ρ, ρ)
    nlsolve(df, ρ; iterations=iterations, ftol=ftol, xtol=xtol, method=method, m=m, kwargs...)
    # nlsolve has a problem: it cannot deal with G0 being a matrix.
    # That's ok for Hartree meanfield, but sucks hard for Fock meanfield.
end

################################################################################
################################################################################

# function solve_selfconsistent(hamiltonian::Function, v::Function,
#     G0::AbstractArray{Float64,N}, ks::AbstractMatrix{Float64}, filling::Float64;
#     mode=:nothing, kwargs...) where N
#
#     if mode==:hartree
#         # Note: this option assumes that v is the Fourier transform of the interaction potential
#         mf_hartree = get_hartree_operator(v)
#         solve_op_selfconsistent(hamiltonian, mf_hartree, G0, ks, filling; kwargs...)
#     else
#         # by default we assume that v is the a single-particle term that implicitly depends on the eigensolution
#         solve_op_selfconsistent(hamiltonian, v, G0, ks, filling; kwargs...)
#     end
# end

################################################################################
################################################################################

function initialize_density(N::Int; mode=:random)
    @assert mod(N,2)==0 # make sure we actually have spin d.o.f.

    if mode==:ferro
        Nhalf = Int(N/2)
        return 2 .* [ones(Float64, Nhalf); zeros(Float64, Nhalf)]
    elseif mode==:antiferro
        Nhalf = Int(N/2)
        return [ones(Float64, Nhalf); -ones(Float64, Nhalf)]
    elseif mode==:zeros
        return zeros(Float64, N)
    elseif mode==:random
        return rand(Float64, N)
    end
end

################################################################################
################################################################################

function get_hartree_operator(v; format=:sparse)
    """
        Builds the mean-field term in the Hamiltonian given the Fourier transform
        of an interaction potential (as e.g. provided by build_meanfield_op.jl).
    """

    function mf_op(n,k)
        vsym(q) = 0.5 .* (v(q).+(v(q))')
        mat = spdiagm(0 => vsym(zero(k)) * n)

        if format==:dense
            mat = Matrix(mat)
        end

        mat
    end
    # function mf_op(ρ,k)
    #     vsym(q) = 0.5 .* (v(q).+(v(q))')
    #
    #     vsym(zero(k)) * ρ
    # end

    return mf_op
end

function build_Hartree(v, args...; kwargs...)
    @warn "build_Hartree(...) has been renamed to get_hartree_operator(...) and will be removed in the future."
    get_hartree_operator(v, args...; kwargs...)
end



function E_Hartree(v0, G0)
    """
        Hartree energy E0 for mean field G0.
        v0: interaction potential at q=0.
    """

    - 1/2 * G0' * v0 * G0
end
