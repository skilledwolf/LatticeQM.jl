
using NLsolve


function solve_op_selfconsistent(hamiltonian::Function, mf_op::Function, G0::AbstractArray{Float64,N}, ks::AbstractMatrix{Float64}, filling::Float64; iterations = 500, ftol=1e-7, xtol=1e-7, method=:anderson, m=5, kwargs...) where N

    type = issparse(hamiltonian(ks[:,1])) ? :sparse : :dense # Decide if the Hamiltonian is sparse

    function f!(G1, G0)
        # Update meanfield Hamiltonian and chemical potential
        new_hamiltonian(k) = hamiltonian(k) .+ mf_op(G0, k)
        μ = chemical_potential(new_hamiltonian, ks, filling; type=type)

        # Obtain the meanfield of the updated Hamiltonian
        density_parallel!(G1, new_hamiltonian, ks, μ)#; format=format)

        G1 .-= G0
        nothing
    end

    df = OnceDifferentiable(f!, G0, G0)
    nlsolve(df, G0; iterations=iterations, ftol=ftol, xtol=xtol, method=method, m=m, kwargs...)
end

################################################################################
################################################################################

function solve_selfconsistent(hamiltonian::Function, v::Function, G0::AbstractArray{Float64,N}, ks::AbstractMatrix{Float64}, filling::Float64; mode=:nothing, kwargs...) where N
    if mode==:hartree
        # Note: this option assumes that v is the Fourier transform of the interaction potential
        mf_hartree = build_Hartree(v)
        solve_op_selfconsistent(hamiltonian, mf_hartree, G0, ks, filling; kwargs...)
    else
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
