
using Distributed
using SharedArrays

using ..Utils: fermidirac

function densitymatrix!(ρ0::AbstractMatrix, ϵs::AbstractVector, U::AbstractMatrix, μ::Float64; T::Float64=0.01, φk::ComplexF64=1.0+0.0im)

    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        ρ0[:] .+= (fermidirac(ϵ-μ; T=T) .* (ψ * ψ') .* φk)[:]
    end

    nothing
end

using ..TightBinding: dim
function densitymatrix(H, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    d = dim(H, ks)
    ρ0 = zeros(ComplexF64, d, d)
    densitymatrix!(ρ0, H, ks, μ; kwargs...)

    ρ0
end

function densitymatrix!(ρ0::AbstractMatrix, H, ks::AbstractMatrix, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    ρ0[:] .= zero(ρ0)[:]
    L = size(ks)[2]

    spectrumH = spectrum(H; kwargs...)

    ρ0[:] = @distributed (+) for j=1:L
        ρtmp = zero(ρ0)    ## <-- it annoys me that I don't know how to get around this allocation

        ϵs, U = spectrumH(ks[:,j])
        densitymatrix!(ρtmp, ϵs, U, μ; T=T) #; format=format)

        ρtmp
    end

    ρ0[:] .= ρ0[:] ./ L

    nothing
end

###################################################################################################
###################################################################################################
###################################################################################################

function densitymatrix!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    L = size(ks,2)

    energies0_k = convert(SharedArray, zeros(Float64, L))
    spectrumf = spectrum(H; kwargs...)

    for (δL,ρ0)=ρs
        ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    @sync @distributed for i_=1:L
        k = ks[:,i_]
        energies_k, U_k = spectrumf(k) #@time

        for δL=keys(ρs)
            for (ϵ, ψ) in zip(energies_k, eachcol(U_k)) # this loop used to be seperate in ρ_k!(...)
                ρs[δL][:] += (fermidirac(real(ϵ)-μ; T=T) .* transpose(ψ * ψ') .* fourierphase(-k, δL))[:]
            end
        end

        energies0_k[i_] = groundstate_sumk(real(energies_k), μ)
    end

    for δL = keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies0_k)/L # return the groundstate energy
end

###################################################################################################
# Legacy maps
###################################################################################################
export ρ_L
@legacyalias densitymatrix ρ_L
