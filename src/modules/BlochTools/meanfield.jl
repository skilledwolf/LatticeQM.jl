
using NLsolve


function solve_op_selfconsistent(hamiltonian::Function, mf_op::Function, G0::AbstractArray{Float64,N}, ks::AbstractMatrix{Float64}, filling::Float64; iterations = 500, ftol=1e-7, xtol=1e-7, method=:anderson, m=5, kwargs...) where N

    type = issparse(hamiltonian(ks[:,1])) ? :sparse : :dense # Decide if the Hamiltonian is sparse

    function f!(G1, G0)
        # Update meanfield Hamiltonian and chemical potential
        new_hamiltonian(k) = hamiltonian(k) .+ mf_op(G0, k)
        μ = chemical_potential(new_hamiltonian, ks, filling; type=type)

        # Obtain the meanfield of the updated Hamiltonian
        density_parallel!(G1, new_hamiltonian, ks, μ)#; format=format)

        # Update meanfield for next iteration step
        G1 .-= G0
        nothing
    end

    df = OnceDifferentiable(f!, G0, G0)
    nlsolve(df, G0; iterations=iterations, ftol=ftol, xtol=xtol, method=method, m=m, kwargs...)
    # nlsolve has a problem: it cannot deal with G0 being a matrix.
    # That's ok for Hartree meanfield, but sucks hard for Fock meanfield.
end

################################################################################
################################################################################

function solve_selfconsistent(hamiltonian::Function, v::Function, G0::AbstractArray{Float64,N}, ks::AbstractMatrix{Float64}, filling::Float64; mode=:nothing, kwargs...) where N
    if mode==:hartree
        # Note: this option assumes that v is the Fourier transform of the interaction potential
        mf_hartree = build_Hartree(v)
        solve_op_selfconsistent(hamiltonian, mf_hartree, G0, ks, filling; kwargs...)
    else
        # by default we assume that v is the a single-particle term that implicitly depends on the eigensolution
        solve_op_selfconsistent(hamiltonian, v, G0, ks, filling; kwargs...)
    end
end

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

function build_Hartree(v; format=:sparse)
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

    return mf_op
end

function E_Hartree(v0, G0)
    """
        Hartree energy E0 for mean field G0.
        v0: interaction potential at q=0.
    """

    - 1/2 * G0' * v0 * G0
end
