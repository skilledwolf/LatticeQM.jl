using Distributed
using ProgressMeter

function ldos!(n::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix, ωs::AbstractVector; Γ::Real=0.1)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag( abs2.(ψ)./(ω + 1.0im * Γ - ϵ) )
        end
    end
    n[:] ./= size(ωs)
end

function ldos!(n::AbstractVector, H, ks::AbstractMatrix, ωs::AbstractVector; Γ::Real=0.1, kwargs...)
    L = size(ks,2)

    spectrumf = spectrum(H; kwargs...)

    n[:] = @sync @showprogress 1 "Computing LDOS..." @distributed (+) for j=1:L
        n0 = zero(n)
        ϵs, U = spectrumf(ks[:,j])
        ldos!(n0, ϵs, U, ωs; Γ=Γ)
        n0
    end
    n[:] .= -n ./ L ./ π
end


using ..TightBinding: dim

ldos(H, ks, frequency::Real; kwargs...) = ldos(H, ks, [frequency]; kwargs...)
function ldos(H, ks, frequencies::AbstractVector; format=:sparse, kwargs...)
    n = zeros(Float64, dim(H,ks))
    ldos!(n, H, ks, frequencies; format=format, kwargs...)
    n
end
