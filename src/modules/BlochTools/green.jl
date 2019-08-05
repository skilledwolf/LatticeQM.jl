
function G!(G0::AbstractMatrix{ComplexF64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; format=:auto)
    G0[:] .= zero(G0)[:]

    # if format==:auto # Decide if the matrix is dense or sparse
    #     format = issparse(hamiltonian(ks[:,1])) ? :sparse : :dense
    # end

    L = size(ks)[2]

    @inbounds for j=1:L # @todo: this should be paralellized
        G_at_k!(G0, hamiltonian, ks[:,j], μ)#; format=format)
    end

    G0[:] .= G0[:] ./ L

    nothing
end

function G_at_k!(G0::AbstractMatrix{ComplexF64}, hamiltonian::Function, k::AbstractVector{Float64}, μ::Float64)
    ϵs, U = spectrum(hamiltonian; format=:dense)(k)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        if ϵ <= μ
            G0[:] .+= (ψ * ψ')[:]
        end
    end
end

function G(hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)

    G0 = zero(Matrix(hamiltonian(ks[:,1])))
    G!(G0, hamiltonian, ks, μ; kwargs...)

    G0
end

using Distributed
# using SharedArrays

function G_parallel!(G0::AbstractMatrix{ComplexF64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)
    G0[:] .= zero(G0)[:]

    # if format==:auto # Decide if the matrix is dense or sparse
    #     format = issparse(hamiltonian(ks[:,1])) ? :sparse : :dense
    # end

    L = size(ks)[2]

    G[:] = @distributed (+) for j=1:L # @todo: this should be paralellized

        G0 = zero(G0)                                ## <-- it annoys me that I don't know how to get around this allocation
        G_at_k!(G0, hamiltonian, ks[:,j], μ)#; format=format)

        G0 .= G0 ./ L
    end

    nothing
end
