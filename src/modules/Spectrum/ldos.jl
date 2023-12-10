using Distributed
using ProgressMeter

import LatticeQM.Eigen
import LatticeQM.Spectrum

const PROGRESSBAR_LDOS_DEFAULTLABEL = "LDOS"::String

function ldos!(n::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix, ωs::AbstractVector; Γ::Real=0.1)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag( abs2.(ψ)./(ω + 1.0im * Γ - ϵ) )
        end
    end
    n[:] ./= size(ωs)
end

function ldos!(n::AbstractVector, H, ks::AbstractMatrix, ωs::AbstractVector; Γ::Real=0.1, progress_label=PROGRESSBAR_LDOS_DEFAULTLABEL, kwargs...)
    L = size(ks,2)


    n[:] = @sync @showprogress dt=Spectrum.PROGRESSBAR_MINTIME desc=progress_label enabled=Spectrum.PROGRESSBAR_SHOWDEFAULT @distributed (+) for j = 1:L
        n0 = zero(n)
        ϵs, U = Eigen.geteigen(H(ks[:, j]); kwargs...)
        ldos!(n0, ϵs, U, ωs; Γ=Γ)
        n0
    end
    n[:] .= -n  ./ L ./ π
end

ldos(H, ks, frequency::Real; kwargs...) = ldos(H, ks, [frequency]; kwargs...)
function ldos(H, ks, frequencies::AbstractVector; format=:sparse, kwargs...)
    n = zeros(Float64, dim(H,ks))
    ldos!(n, H, ks, frequencies; format=format, kwargs...)
    n
end
