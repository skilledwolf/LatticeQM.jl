
using NLsolve


function update_density!(G1, G0, hamiltonian::Function, mf_op::Function, ks::AbstractMatrix{Float64}, filling::Float64; mode=:sparse)

    # G0 = reshape(G0, (N,N)) # very sub-optimal

    new_hamiltonian(k) = hamiltonian(k) .+ mf_op(G0, k)
    μ = chemical_potential(hamiltonian, ks, filling; mode=mode)

    # G1 = reshape(G1) # very sub-optimal
    G1[:] .= zero(G1)
    density!(G1, new_hamiltonian, ks, μ) #density_parallel!

    nothing
end


function solve_op_selfconsistent(hamiltonian::Function, mf_op::Function, ks::AbstractMatrix{Float64}, G0::AbstractArray{Float64,N}, filling::Float64; iterations = 500, ftol=1e-7, xtol=1e-7, method=:anderson, m=5, kwargs...) where N

    if issparse(hamiltonian(ks[:,1]))
        mode=:sparse
    else
        mode=:dense
    end

    function f!(F,x)
        update_density!(F, x, hamiltonian, mf_op, ks, filling; mode=mode)
        F .-= x
        nothing
    end

    df = OnceDifferentiable(f!, G0, G0)

    nlsolve(df, G0; iterations=iterations, ftol=ftol, xtol=xtol, method=method, m=m, kwargs...)
end

function solve_op_selfconsistent(hamiltonian::Function, mf_op::Function, G0::AbstractArray{Float64,N}, filling::Float64; nk::Int=100, d::Int=2, kwargs...) where N

    ks = rand(Float64, (d,nk)) # random points in the Brillouin zone (which is assumed to be mapped onto a unit cube)

    solve_op_selfconsistent(hamiltonian, mf_op, ks, G0, filling; kwargs...)
end

function solve_selfconsistent(hamiltonian::Function, v::Function, G0::AbstractArray{Float64,N}, filling::Float64; mode=:hartree, nk::Int=100, d::Int=2, kwargs...) where N
    if mode==:hartree
        mf_hartree = build_Hartree(v)
        solve_op_selfconsistent(hamiltonian, mf_hartree, G0, filling; nk=nk, d=d, kwargs...)
    else
        error("Requested mean-field operator not implemented.")
    end
end

# using SharedArrays
#
# function solve_selfconsistent_parallel(hamiltonian::Function, v::Function, G0::AbstractArray{Float64,N}, filling::Float64; mode=:hartree, nk::Int=100, d::Int=2, kwargs...) where N
#     if mode==:hartree
#         mf_hartree = build_Hartree(v)
#         G0 = convert(SharedArray, G0)
#         solve_op_selfconsistent(hamiltonian, mf_hartree, G0, filling; nk=nk, d=d, kwargs...)
#     else
#         error("Requested mean-field operator not implemented.")
#     end
# end


function build_Hartree(v)

    function mf_op(n,k)
        vsym(q) = 0.5 .* (v(q).+(v(q))')
        spdiagm(0 => vsym(zero(k)) * n)
    end

    return mf_op
end


function initialize_density(N::Int; mode=:random)
    @assert mod(N,2)==0 # make sure we actually have spin d.o.f.

    if mode==:ferro
        Nhalf = Int(N/2)
        return 2 .* [ones(Float64, Nhalf) zeros(Float64, Nhalf)]
    elseif mode==:zeros
        return zeros(Float64, N)
    elseif mode==:random
        return random(Float64, N)
    end
end
