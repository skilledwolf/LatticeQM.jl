using Distributed

function ldos_at_k!(n::AbstractVector{Float64}, eigen::Function, k::AbstractVector{Float64}, ωs::AbstractVector{Float64}; Γ::Float64=0.1)
    ϵs, U = eigen(k)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag( abs2.(ψ)./(ω + 1.0im * Γ - ϵ) )
        end
    end

    nothing
end

function ldos_parallel!(n::AbstractVector{Float64}, eigen::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64=0.1)

    L = size(ks)[2]

    n[:] = @distributed (+) for j=1:L

        n0 = zero(n)                                ## <-- it annoys me that I don't know how to get around this allocation
        ldos_at_k!(n0, eigen, ks[:,j], ωs; Γ=Γ)

        n0
    end

    n[:] .= -n[:] ./ L ./ π

    nothing
end

function ldos(hamiltonian::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64, type=:sparse, kwargs...)

    n = zeros(Float64, size(hamiltonian(ks[:,1]))[1])

    if type==:sparse
        eigen = eigen_sparse(hamiltonian; kwargs...)
    elseif type==:dense
        eigen = eigen_dense(hamiltonian; kwargs...)
    else
        error("Invalid type-mode in ldos.")
    end

    ldos_parallel!(n, eigen, ks, ωs; Γ=Γ)

    n
end
