
using Distributed
using SharedArrays
using ProgressMeter

using ..Utils: fermidirac

using ..TightBinding: dim
function densitymatrix(H, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    d = dim(H, ks)
    ρ0 = zeros(ComplexF64, d, d)
    densitymatrix!(ρ0, H, ks, μ; kwargs...)

    ρ0
end

# function densitymatrix!(ρ0::AbstractMatrix, H, ks::AbstractMatrix, μ::Float64=0.0; T::Float64=0.01, kwargs...)
#     ρ0[:] .= zero(ρ0)[:]
#     L = size(ks)[2]
#
#     spectrumH = spectrum(H; kwargs...)
#
#     ρ0[:] = @distributed (+) for j=1:L
#         ρtmp = zero(ρ0)    ## <-- it annoys me that I don't know how to get around this allocation
#
#         ϵs, U = spectrumH(ks[:,j])
#         densitymatrix!(ρtmp, ϵs, U, μ; T=T) #; format=format)
#
#         ρtmp
#     end
#
#     ρ0[:] .= ρ0[:] ./ L
#
#     nothing
# end
densitymatrix(ϵ::Number, ψ::AbstractVector; T::Float64=0.01) = fermidirac(real(ϵ); T=T) .* transpose(ψ * ψ') #transpose(ψ * ψ') # (ψ * ψ')

function densitymatrix!(ρ0::AbstractMatrix, ϵs::AbstractVector, U::AbstractMatrix; φk::ComplexF64=1.0+0.0im, kwargs...)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        # ρ0[:,:] .+= (densitymatrix(ϵ, ψ; kwargs...) .* φk)
        broadcast!(+, ρ0, ρ0, (densitymatrix(ϵ, ψ; kwargs...) .* φk))
    end
    ρ0
end

function densitymatrix!(ρ0::AbstractMatrix, δL::AbstractVector, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; kwargs...)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        # ρ0[:,:] .+= (densitymatrix(ϵ, ψ; kwargs...) .* fourierphase(-k, δL)) # ϵ-μ # -k
        broadcast!(+, ρ0, ρ0, (densitymatrix(ϵ, ψ; kwargs...) .* fourierphase(-k, δL)))
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

# function densitymatrix_parallel!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
#     L = size(ks,2)

#     energies0_k = zeros(Float64, L) #convert(SharedArray, zeros(Float64, L))
#     spectrumf = spectrum(H; kwargs...)

#     for (δL,ρ0)=ρs
#         ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
#     end

#     channel = RemoteChannel(()->Channel{Tuple{Int, Vector{Float64}, Matrix{Complex}}}(L), 1)
#     @sync begin
#         @async begin # update ρs
#             done = 0
#             while done < L
#                 (i_, ϵs, U) = take!(channel) # read the result from channel (wait if necessary)
#                 densitymatrix!(ρs, ks[:,i_], ϵs.-μ, U; T=T)
#                 energies0_k[i_] = groundstate_sumk(real(ϵs), μ)
#                 done = done+1
#             end
#         end

#         @async begin # compute spectrum at different k points asynchronosly (good for large/huge systems)
#             # @sync @showprogress 1 "Eigensolver... " @distributed for i_=1:L
#             @sync @showprogress 10 "Eigensolver... " @distributed for i_=1:L
#                 k = ks[:,i_]
#                 energies_k, U_k = spectrumf(k) # calculation
#                 put!(channel, (i_, real.(energies_k), U_k)) # passing the result to the channel
#             end
#         end
#     end

#     for δL = keys(ρs)
#         ρs[δL][:] ./= L
#     end

#     sum(energies0_k)/L # return the groundstate energy
# end

using ..TightBinding: efficientformat, efficientzero, flexibleformat!
using ProgressBars

function densitymatrix_multithread!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    L = size(ks,2)

    energies = zeros(Float64, L)
    spectrumf = spectrum(H; kwargs...)

    for (δL,ρ0)=ρs
        ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    Threads.@threads for i_=ProgressBar(1:L)
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        lock(ρs) do
            for δL=keys(ρs)
                densitymatrix!(ρs[δL], δL, k, ϵs.-μ, U; T=T)
            end
        end

        # densitymatrix!(ρs, k, ϵs.-μ, U; T=T)
        energies[i_] = groundstate_sumk(real(ϵs), μ)
    end

    for δL = keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies)/L # return the groundstate energy
end

function densitymatrix_parallel!(ρs::Dict{Vector{Int},SharedMatrix}, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    L = size(ks,2)

    energies = SharedArray(zeros(Float64, L))

    for δL=keys(ρs)
        ρs[δL][:] .= 0
    end

    spectrumf = spectrum(H; kwargs...)

    @sync @showprogress 10 "Eigensolver... " @distributed for i_=1:L
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        for δL=keys(ρs)
            densitymatrix!(view(ρs[δL],:,:), δL, k, ϵs.-μ, U; T=T)
        end

        energies[i_] = groundstate_sumk(real(ϵs), μ)
    end

    for δL=keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies)/L # return the groundstate energy
end

# function densitymatrix_parallel!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
#     L = size(ks,2)

#     energies = SharedArray(zeros(Float64, L))
#     spectrumf = spectrum(H; kwargs...)

#     zeromat, δLs = efficientzero(ρs)

#     ρsMat = @sync @showprogress 10 "Eigensolver... " @distributed (+) for i_=1:L
#         k = ks[:,i_]
#         ϵs, U = spectrumf(k) #@time

#         ρ0 = deepcopy(zeromat)
#         for (j_,δL)=enumerate(δLs)
#             densitymatrix!(view(ρ0,:,:,j_), δL, k, ϵs.-μ, U; T=T)
#             # densitymatrix!(view(ρsMat, :,:,j_), δL, k, ϵs.-μ, U; T=T)
#         end

#         energies[i_] = groundstate_sumk(real(ϵs), μ)

#         ρ0
#     end

#     flexibleformat!(ρs, ρsMat/L, δLs)

#     sum(energies)/L # return the groundstate energy
# end

function densitymatrix_serial!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01, kwargs...)
    L = size(ks,2)

    energies = zeros(Float64, L)
    spectrumf = spectrum(H; kwargs...)

    for (δL,ρ0)=ρs
        ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    @showprogress 10 "Eigensolver... " for i_=1:L
        k = ks[:,i_]
        ϵs, U = spectrumf(k) #@time

        for δL=keys(ρs)
            densitymatrix!(ρs[δL], δL, k, ϵs.-μ, U; T=T)
        end

        # densitymatrix!(ρs, k, ϵs.-μ, U; T=T)
        energies[i_] = groundstate_sumk(real(ϵs), μ)
    end

    for δL = keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies)/L # return the groundstate energy
end


function densitymatrix!(ρs::AnyHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; multimode=:serial, kwargs...)

    if multimode==:parallel && nprocs()>1

        densitymatrix_parallel!(ρs, H, ks, μ; kwargs...)
    elseif multimode==:multithread && Threads.nthreads()>1
        densitymatrix_multithread!(ρs, H, ks, μ; kwargs...)
    else
        densitymatrix_serial!(ρs, H, ks, μ; kwargs...)
    end
end

