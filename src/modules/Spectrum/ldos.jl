using Distributed

function ldos_at_k!(n::AbstractVector{Float64}, spectrum::Function, k::AbstractVector{Float64}, ωs::AbstractVector{Float64}; Γ::Float64=0.1)
    ϵs, U = spectrum(k)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag( abs2.(ψ)./(ω + 1.0im * Γ - ϵ) )
        end
    end
    n[:] ./= size(ωs)

    nothing
end

function ldos!(n::AbstractVector{Float64}, spectrum::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64=0.1)

    L = size(ks,2)

    n[:] = @distributed (+) for j=1:L

        n0 = zero(n)                                ## <-- it annoys me that I don't know how to get around this allocation
        ldos_at_k!(n0, spectrum, ks[:,j], ωs; Γ=Γ)

        n0
    end

    n[:] .= -n ./ L ./ π

    nothing
end

function ldos(hamiltonian::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64, type=:sparse, kwargs...)

    n = zeros(Float64, size(hamiltonian(ks[:,1]))[1])

    eigen = spectrum(hamiltonian; type=type, kwargs...)

    ldos!(n, eigen, ks, ωs; Γ=Γ)

    n
end

ldos(hamiltonian::Function, ks::AbstractMatrix{Float64}, ω::Float64; kwargs...) = ldos(hamiltonian, ks, [ω]; kwargs...)
