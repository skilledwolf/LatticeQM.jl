using Distributed

using ProgressMeter

function ldos!(n::AbstractVector{Float64}, spectrum::Function, k::AbstractVector{Float64}, ωs::AbstractVector{Float64}; Γ::Float64=0.1)
    ϵs, U = spectrum(k)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag( abs2.(ψ)./(ω + 1.0im * Γ - ϵ) )
        end
    end
    n[:] ./= size(ωs)
end

function ldos!(n::AbstractVector{Float64}, spectrum::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64=0.1)
    L = size(ks,2)

    n[:] = @sync @showprogress 1 "Computing LDOS..." @distributed (+) for j=1:L
        n0 = zero(n)
        ldos!(n0, spectrum, ks[:,j], ωs; Γ=Γ)
        n0
    end
    n[:] .= -n ./ L ./ π
end

ldos(hamiltonian::Function, ks::AbstractMatrix{Float64}, ω::Float64; kwargs...) = ldos(hamiltonian, ks, [ω]; kwargs...)
function ldos(hamiltonian::Function, ks::AbstractMatrix{Float64}, ωs::AbstractVector{Float64}; Γ::Float64, format=:sparse, kwargs...)

    n = zeros(Float64, size(hamiltonian(ks[:,1]))[1])

    eigen = spectrum(hamiltonian; format=format, kwargs...)

    ldos!(n, eigen, ks, ωs; Γ=Γ)
    n
end
