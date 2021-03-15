@fastmath function groundstate_sumk(ϵs_k::AbstractVector{T}, μ::T=0.0) where {T<:Number}
    tmp = 0.0
    for ϵ in ϵs_k
        if ϵ <= μ
            tmp += ϵ
        end
    end

    tmp
end

function groundstate_energy(ϵs::Function, ks::AbstractMatrix{T}, μ::T=0.0; kwargs...) where {T<:Number}
    # Σ = ϵs(hamiltonian; format=format)
    L = size(ks)[2]

    ϵGS = @distributed (+) for j=1:L
        groundstate_sumk(ϵs(ks[:,j]), μ)
    end

    ϵGS / L
end

# function groundstate_energy_multithread(ϵs::Function, ks::AbstractMatrix{T}, μ::T=0.0; kwargs...) where {T<:Number}
#     # Σ = ϵs(hamiltonian; format=format)
#     L = size(ks)[2]

#     ϵGS = Threads.Atomic{Float64}(0)
#     Threads.@threads for j=1:L
#         e0 = groundstate_sumk(ϵs(ks[:,j]), μ)
#         Threads.atomic_add!(ϵGS, e0)
#     end

#     ϵGS / L
# end
