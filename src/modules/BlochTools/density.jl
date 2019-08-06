using Distributed

function density_at_k!(n::AbstractVector{Float64}, spectrum_k, μ::Float64)
    ϵs, U = spectrum_k
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        if ϵ <= μ
            n[:] .+= abs2.(ψ)
        end
    end
end

function density(hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; format=:dense, kwargs...)
    Σ = spectrum(hamiltonian; format=format)

    n = zeros(Float64, size(hamiltonian(ks[:,1]))[1])
    density!(n, Σ, ks, μ; kwargs...)

    n
end

function density!(n::AbstractVector{Float64}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)
    n[:] .= zero(n)

    L = size(ks)[2]

    n[:] = @distributed (+) for j=1:L # @todo: this should be paralellized

        n0 = zero(n)    ## <-- it annoys me that I don't know how to get around this allocation
        density_at_k!(n0, spectrum(ks[:,j]), μ)

        n0 .= n0 ./ L
    end

    nothing
end

# function density!(n::AbstractVector{Float64}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; format=:auto)
#     n[:] .= zero(n)
#
#     # if format==:auto # Decide if the matrix is dense or sparse
#     #     format = issparse(hamiltonian(ks[:,1])) ? :sparse : :dense
#     # end
#
#     L = size(ks)[2]
#
#     @inbounds for j=1:L # @todo: this should be paralellized
#         density_at_k!(n, spectrum(ks[:,j]), μ)#; format=format)
#     end
#
#     n[:] .= n[:] ./ L
#
#     nothing
# end
