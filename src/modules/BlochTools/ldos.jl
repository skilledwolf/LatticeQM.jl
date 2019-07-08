using Distributed

function ldos_at_k!(n::AbstractVector{Float64}, hamiltonian::Function, k::AbstractVector{Float64}, ωs::AbstractVector{Float64}; Γ::Float64)
    ϵs, U = eigen_dense(hamiltonian)(k)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag( abs2.(ψ)./(ω + 1.0im * Γ - ϵ) )
        end
    end

    nothing
end


function ldos(hamiltonian::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64, kwargs...)

    n = zeros(Float64, size(hamiltonian(ks[:,1]))[1])
    ldos_parallel!(n, hamiltonian, ks, ωs::AbstractVector{Float64}; Γ::Float64, kwargs...)

    n
end


function ldos_parallel!(n::AbstractVector{Float64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64)

    L = size(ks)[2]

    n[:] = @distributed (+) for j=1:L # @todo: this should be paralellized

        n0 = zero(n)                                ## <-- it annoys me that I don't know how to get around this allocation
        ldos_at_k!(n0, hamiltonian, ks[:,j], ωs, Γ)

        n0
    end

    n[:] .= -n[:] ./ L ./ π

    nothing
end
