
using Distributed
using SharedArrays
using ProgressMeter

import ..Utils: fermidirac
import ..TightBinding: Hops, AnyHops, dim

function getdensitymatrix(H, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    d = dim(H, ks)
    ρ0 = zeros(ComplexF64, d, d)
    densitymatrix!(ρ0, H, ks, μ; kwargs...)

    ρ0
end

densitymatrix(ϵ::Number, ψ::AbstractVector; T::Real=0.01) = fermidirac(real(ϵ); T=T) .* transpose(ψ * ψ') #transpose(ψ * ψ') # (ψ * ψ')

function densitymatrix!(ρ0::AbstractMatrix, ϵs::AbstractVector, U::AbstractMatrix; φk::ComplexF64=1.0+0.0im, T=0.01, kwargs...)
    fd = fermidirac.(real.(ϵs); T=T)

    for m in 1:length(ϵs)
        ρ0[:,:] .+= (fd[m] .* φk .* (conj.(U[:,m]) * transpose(U[:,m])))
    end
    ρ0
end


import ..TightBinding: fourierphase

function densitymatrix!(ρ0::AbstractMatrix, δL::AbstractVector, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; T=0.01, kwargs...)
    
    phase = fourierphase(-k, δL)
    densitymatrix!(ρ0, ϵs, U; φk=phase, T=T, kwargs...)
end

function mapdensitymatrix!(ρ0::AbstractMatrix, δLs::AbstractVector{AbstractVector}, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; T=0.01, kwargs...)

    M = zero(ρ0)
    densitymatrix!(M, ϵs, U; φk=phase, T=T, kwargs...)

    for δL=δLs 
        ρ0[:,:] .+= M .* fourierphase(-k, δL)
    end

    ρ0
end

function densitymatrix!(ρs::AnyHops, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; kwargs...)
 
    for δL=keys(ρs)
        densitymatrix!(ρs[δL], δL, k, ϵs, U; kwargs...)
    end
    ρs
end

###################################################################################################
###################################################################################################
###################################################################################################
using ProgressBars

import ..Spectrum: spectrum, groundstate_sumk
import ..TightBinding: efficientzero, flexibleformat!, fourierphase

function densitymatrix_multithread!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks,2)

    energies = zeros(Float64, L)
    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    for (δL,ρ0)=ρs
        ρs[δL][:] .= 0 #convert(SharedArray, zero(ρ0))[:]
    end

    Msize = size(first(values(ρs)))
    Mtype = eltype(first(values(ρs)))

    lk = Threads.ReentrantLock()
    Threads.@threads for i_=ProgressBar(1:L)#i_=1:L
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        M = zeros(Mtype, Msize)
        densitymatrix!(M, ϵs.-μ, U; T=T)

        lock(lk) do 
            for δL=keys(ρs)
                    ρs[δL][:,:] .+= (M .* fourierphase(k, δL))
            end
        end

        energies[i_] = groundstate_sumk(real(ϵs), μ)
    end

    for δL = keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies)/L # return the groundstate energy
end


using ..TightBinding: Hops, AnyHops

import ..Spectrum: spectrum, groundstate_sumk
import ..TightBinding: efficientzero, flexibleformat!, fourierphase

#### Nice idea but less performant... large overhead for copying data between processes
# function densitymatrix_pmap_async!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
#     L = size(ks,2)

#     energies = zeros(Float64, L)

#     function spectrumf(k)
#         spectrum(H(k); kwargs...)
#     end

#     for (δL,ρ0)=ρs
#         ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
#     end

#     zeromat = complex(zeros(size(first(values(ρs)))))

#     channel = RemoteChannel(()->Channel{Tuple{Int, Vector{Float64}, Matrix{Complex}}}(L))
#     @sync begin

#         @async begin # update G
#             done = 0
#             wait(channel)
#             while done < L
#                 (i_, ϵs, M) = take!(channel) # read the result from channel (wait if necessary)

#                 for δL=keys(ρs)
#                     ρs[δL][:,:] .+= (M .* fourierphase(ks[:,i_], δL))
#                 end

#                 energies[i_] = groundstate_sumk(real(ϵs), μ)
#                 done = done+1
#             end
#         end

#         @async begin # compute spectrum at different k points asynchronosly (good for large/huge systems)
#             pmap(1:L) do i_
#                 k = ks[:,i_]
#                 # println("Calculating spectrum: $i_")
#                 ϵs, U = spectrumf(k) # calculation

#                 # println("Calculating spectrum: $i_")
#                 M = zero(zeromat)
#                 densitymatrix!(M, ϵs.-μ, U; T=T)

#                 # println("Sending to channel: $i_")
#                 put!(channel, (i_, real.(ϵs), M)) # passing the result to the channel
#                 nothing
#             end
#         end
#     end

#     for δL = keys(ρs)
#         ρs[δL][:] ./= L
#     end

#     sum(real(energies))/L # return the groundstate energy
# end

# function densitymatrix_pmap!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
#     L = size(ks,2)

#     function spectrumf(k)
#         spectrum(H(k); kwargs...)
#     end

#     Ls = keys(ρs)
#     M0 = first(values(ρs))

#     res = pmap(eachcol(ks)) do k
#         ϵs, U = spectrumf(k) 

#         M = zero(M0)
#         densitymatrix!(M, ϵs.-μ, U; T=T)

#         ρ0 = Dict(L => M .* fourierphase(-k, δL) for L in Ls)

#         ρ0, groundstate_sumk(real(ϵs), μ)
#     end

#     ρsnew = mergewith(+, (x[1] for x=res)...)
#     energies = sum(x[2] for x=res)/L

#     for L in Ls
#         ρs[L] .= ρsnew[L]
#     end

#     energies # return the groundstate energy
# end

function densitymatrix_pmap!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks,2)

    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    zeromat, δLs = efficientzero(ρs)

    ρ0 = SharedArray(zeromat)

    energy = sum(pmap(1:L) do i_
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        M = zero(ρ0[:,:,1])

        densitymatrix!(M, ϵs.-μ, U; T=T)

        for (j_,δL)=enumerate(δLs)
            ρ0[:,:,j_] .+= (M .* fourierphase(-k, δL))
        end

        groundstate_sumk(real(ϵs), μ)
    end)/L

    flexibleformat!(ρs, ρ0/L, δLs)

    energy # return the groundstate energy
end

# function densitymatrix_pmap!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
#     L = size(ks,2)

#     energies = SharedArray(zeros(Float64, L))
#     function spectrumf(k)
#         spectrum(H(k); kwargs...)
#     end

#     zeromat, δLs = efficientzero(ρs)

#     ρsMat = sum(pmap(1:L) do i_
#         k = ks[:,i_]
#         ϵs, U = spectrumf(k) #@time

#         ρ0 = zero(zeromat)
#         M = zero(ρ0[:,:,1])

#         densitymatrix!(M, ϵs.-μ, U; T=T)

#         for (j_,δL)=enumerate(δLs)
#             ρ0[:,:,j_] .+= (M .* fourierphase(-k, δL))
#         end

#         energies[i_] = groundstate_sumk(real(ϵs), μ)

#         ρ0
#     end)

#     flexibleformat!(ρs, ρsMat/L, δLs)

#     sum(energies)/L # return the groundstate energy
# end

# deprecated in favor of densitymatrix_pmap
function densitymatrix_distributed!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks,2)

    energies = SharedArray(zeros(Float64, L))
    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    zeromat, δLs = efficientzero(ρs)

    ρsMat = @sync @showprogress progressmin "Eigensolver... " @distributed (+) for i_=1:L
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        ρ0 = zero(zeromat)
        M = zero(ρ0[:,:,1])

        densitymatrix!(M, ϵs.-μ, U; T=T)

        for (j_,δL)=enumerate(δLs)
            ρ0[:,:,j_] .+= (M .* fourierphase(-k, δL))
        end

        energies[i_] = groundstate_sumk(real(ϵs), μ)

        ρ0
    end

    flexibleformat!(ρs, ρsMat/L, δLs)

    sum(energies)/L # return the groundstate energy
end

function densitymatrix_serial!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks,2)

    energies = zeros(Float64, L)
    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    for (δL,ρ0)=ρs
        ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    @showprogress progressmin "Eigensolver... " for i_=1:L
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        M = zero(first(values(ρs)))
        densitymatrix!(M, ϵs.-μ, U; T=T)

        for δL=keys(ρs)
            ρs[δL][:,:] .+= (M .* fourierphase(k, δL))
        end

        energies[i_] = groundstate_sumk(real(ϵs), μ)
    end

    for δL = keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies)/L # return the groundstate energy
end


function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; multimode=:serial, kwargs...)

    if multimode==:distributed && nprocs()>1
        densitymatrix_pmap!(ρs, H, ks, μ; kwargs...)
    elseif multimode==:multithread && Threads.nthreads()>1
        densitymatrix_multithread!(ρs, H, ks, μ; kwargs...)
    else
        densitymatrix_serial!(ρs, H, ks, μ; kwargs...)
    end
end

