#####################################
#
# The implementation in this file needs reviewing and testing!
#
#####################################

using Distributed
using SharedArrays
using ProgressMeter

import ..Utils: fermidirac
using ..TightBinding: Hops, AbstractHops, dim

function green(H, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    d = dim(H, ks)
    G = zeros(ComplexF64, d, d)
    green!(G, H, ks, μ; kwargs...)
    G
end

green(ϵ::Number, ψ::AbstractVector, ω::Number=0.0) = transpose(ψ * ψ') / (ω-ϵ) #transpose(ψ * ψ') # (ψ * ψ')

function green!(G::AbstractMatrix, ϵs::AbstractVector, U::AbstractMatrix; φk::ComplexF64=1.0+0.0im, kwargs...)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        G[:] .+= (green(ϵ, ψ; kwargs...) .* φk)[:]
    end
    G
end

function green!(G::AbstractHops, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; kwargs...)
    for δL=keys(G)
        for (ϵ, ψ) in zip(ϵs, eachcol(U))
            G[δL][:] .+= (green(ϵ, ψ; kwargs...) .* fourierphase(-k, δL))[:] # ϵ-μ # -k
        end
    end
    G
end

###################################################################################################
###################################################################################################
###################################################################################################

import ..Spectrum: spectrum, groundstate_sumk

function green_parallel!(G::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    L = size(ks,2)

    energies0_k = zeros(Float64, L) #convert(SharedArray, zeros(Float64, L))
    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    for δL = keys(G)
        G[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    channel = RemoteChannel(()->Channel{Tuple{Int, Vector{Float64}, Matrix{Complex}}}(L), 1)
    @sync begin
        @async begin # update G
            done = 0
            while done < L
                (i_, ϵs, U) = take!(channel) # read the result from channel (wait if necessary)
                green!(G, ks[:,i_], ϵs.-μ, U; T=T)
                energies0_k[i_] = groundstate_sumk(real(ϵs), μ)
                done = done+1
            end
        end

        @async begin # compute spectrum at different k points asynchronosly (good for large/huge systems)
            @sync @showprogress 1 "Eigensolver... " @distributed for i_=1:L
                k = ks[:,i_]
                energies_k, U_k = spectrumf(k) # calculation
                put!(channel, (i_, real.(energies_k), U_k)) # passing the result to the channel
            end
        end
    end

    for δL = keys(G)
        G[δL][:] ./= L
    end

    sum(energies0_k)/L # return the groundstate energy
end

function green_serial!(G::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    L = size(ks,2)

    energies0_k = zeros(Float64, L)
    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    for δL = keys(G)
        G[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    @showprogress 1 "Eigensolver... " for i_=1:L
        k = ks[:,i_]
        energies_k, U_k = spectrumf(k) #@time

        green!(G, k, energies_k.-μ, U_k; T=T)
        energies0_k[i_] = groundstate_sumk(real(energies_k), μ)
    end

    for δL = keys(G)
        G[δL][:] ./= L
    end

    sum(energies0_k)/L # return the groundstate energy
end

green!(G::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; parallel::Bool=false, kwargs...) = (parallel && nprocs()>1) ? green_parallel!(G, H, ks, μ; kwargs...) : green_serial!(G, H, ks, μ; kwargs...)

