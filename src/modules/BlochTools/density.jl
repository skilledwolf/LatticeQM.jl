
function density!(n::AbstractVector{Float64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)

    L = size(ks)[2]

    @inbounds @simd for j=1:L # @todo: this should be paralellized
        density_at_k!(n, hamiltonian, ks[:,j], μ)
    end

    nothing
end

function density_at_k!(n::AbstractVector{Float64}, hamiltonian::Function, k::AbstractVector{Float64}, μ::Float64)
    ϵs, U = eigen_dense(hamiltonian)(k)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        if ϵ <= μ
            n[:] .+= abs2.(ψ)
        end
    end
end

function density(hamiltonian::Function, ks::AbstractMatrix{Float64}; kwargs...)

    n = zeros(Float64, size(hamiltonian(ks[:,1]))[1])
    density!(n, hamiltonian, ks; kwargs...)

    n
end

# using Distributed
# using SharedArrays

function density_parallel!(n::AbstractVector{Float64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)

    L = size(ks)[2]

    n[:] = @distributed (+) for j=1:L # @todo: this should be paralellized
        ϵs, U = eigen_dense(hamiltonian)(ks[:,j])

        n0 = zero(n)                        ## <-- it annoys me that I don't know how to get around this allocation
        for (ϵ, ψ) in zip(ϵs, eachcol(U))
            if ϵ <= μ
                n0 .+= abs2.(ψ)
            end
        end

        n0
    end

    nothing
end

# function density_parallel!(n::AbstractVector{Float64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)
#
#     L = size(ks)[2]
#
#     Threads.@threads for j=1:L # @todo: this should be paralellized
#         ϵs, U = eigen_dense(hamiltonian)(ks[:,j])
#
#         lock()
#         for (ϵ, ψ) in zip(ϵs, eachcol(U))
#             if ϵ <= μ
#
#                 n[:] .+= abs2.(ψ)
#
#             end
#         end
#         unlock()
#
#     end
#
#     nothing
# end

# using Distributed
# using SharedArray
#
# function density!(n::AbstractVector{Float64}, hamiltonian::Function, ks::AbstractMatrix{Float64}; μ::Float64=0.0)
#
#     ϵ_ψ = energies_wfs(hamiltonian, ks) # get the iterator
#
#     n = convert(SharedArray, n)
#
#     @sync @distributed for F in ϵ_ψ
#         for (ϵ, ψ) in zip(F.values, eachcol(F.vectors))
#             if ϵ < μ
#                 n[:] .+= abs2.(ψ)
#             end
#         end
#     end
#
#     n = convert(Array, n)
#
#     nothing
# end

# function density_filling!(n::AbstractVector{Float64}, hamiltonian::Function, ks::AbstractMatrix{Float64}; filling::Float64)
#
#     μ = chemical_potential(hamiltonian, ks, filling)
#     density_filling!(n, hamiltonian, ks; μ=μ)
#
#     ϵ_ψ = energies_wfs(hamiltonian, ks) # get the iterator
#
#     for F in ϵ_ψ # this loop can (and should) be parallelized
#         for (ϵ, ψ) in zip(F.values, eachcol(F.vectors))
#             if ϵ < μ
#                 n[:] .+= abs2.(ψ)
#             end
#         end
#     end
#
#     nothing
# end
