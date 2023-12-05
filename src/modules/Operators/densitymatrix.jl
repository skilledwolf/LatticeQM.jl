using Distributed
using SharedArrays
using ProgressMeter
import Tullio # for tensor contractions
# import TensorOperations # for tensor contractions

import ..Utils: fermidirac
import ..TightBinding: Hops, AbstractHops, dim
import ..TightBinding: fourierphase
import ..Structure: Mesh, meshweights

densitymatrix(ϵ::Number, ψ::AbstractVector; T::Real = 0.01) = fermidirac(real(ϵ); T = T) .* transpose(ψ * ψ') #transpose(ψ * ψ') # (ψ * ψ')

function densitymatrix!(ρ0::AbstractMatrix, ϵs::AbstractVector, U::AbstractMatrix; φk::ComplexF64 = 1.0 + 0.0im, T = 0.01, kwargs...)
    fd = fermidirac.(real.(ϵs); T = T)

    # note Oct 20 2021: (conj.(U[:,m]) * transpose(U[:,m]))) --> U[:,m] * U[:,m]'))
    Tullio.@tullio ρ0[i, j] += fd[m] * φk * U[i, m] * conj(U[j, m])
    # TensorOperations.@tensoropt ρ0[i, j] += fd[m] * φk * U[i, m] * conj(U[j, m])
    ρ0
end

function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; kwargs...)
    L = size(ks, 2); kweights = fill(1/L, L)
    getdensitymatrix!(ρs, H, ks, kweights, μ; kwargs...)
end

function getdensitymatrix!(ρs::Hops, H, ks::Mesh, μ::Float64=0.0; kwargs...)
    ks = kgrid.points
    kweights = meshweights(kgrid)
    getdensitymatrix!(ρs, H, ks, kweights, μ; kwargs...)
end

function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; multimode=:serial, kwargs...)

    if multimode == :distributed && nprocs() > 1
        densitymatrix_distributed!(ρs, H, ks, kweights, μ; kwargs...) 
        # densitymatrix_dagger!(ρs, H, ks, kweights, μ; kwargs...) 
    else
        densitymatrix_serial!(ρs, H, ks, kweights, μ; kwargs...)
    end
end

###################################################################################################
###################################################################################################
###################################################################################################
# using ProgressBars

import ..Spectrum: spectrum
import ..TightBinding: efficientzero, flexibleformat!, fourierphase

function densitymatrix_distributed!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks, 2)

    energies = SharedArray(zeros(Float64, L))
    δLs = collect(keys(ρs))

    M = length(δLs); Ns = size(H)
    ρsMat = SharedArray(zeros(ComplexF64, Ns..., M))

    @sync @distributed for i_ = 1:L
        k = ks[:, i_]
        ϵs, U = spectrum(H(k); kwargs...)
        fd = fermidirac.(real.(ϵs .- μ); T=T)

        w = kweights[i_]
        phases = [fourierphase(-k, δL) for δL in δLs]

        fdw = fd .* w

        Tullio.@tullio ρsMat[i, j, n] += fdw[m] * phases[n] * U[i, m] * conj(U[j, m])
        # ρsMat .= ρsMat .+ TensorOperations.@tensoropt mat[i, j, n] := fdw[m] * phases[n] * U[i, m] * conj(U[j, m])
        U = nothing

        energies[i_] = real(w * sum(ϵs .* fermidirac(ϵs .- μ)))
    end

    flexibleformat!(ρs, ρsMat, δLs)
    sum(energies) # return the kinetic part of the gs energy
end

# function densitymatrix_distributed!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
#     L = size(ks, 2)

#     energies = SharedArray(zeros(Float64, L))
#     δLs = collect(keys(ρs))

#     ρsMat = @distributed (+) for i_ = 1:L
#         k = ks[:, i_]
#         ϵs, U = spectrum(H(k); kwargs...)
#         fd = fermidirac.(real.(ϵs .- μ); T=T)

#         w = kweights[i_]
#         phases = [fourierphase(-k, δL) for δL in δLs]
        
#         @time Tullio.@tullio cuda = false grad = false ρ0[i, j, n] := fd[m] * phases[n] * U[i, m] * conj(U[j, m])
#         energies[i_] = real(w * sum(ϵs .* fermidirac(ϵs .- μ)))

#         w .* ρ0
#     end

#     flexibleformat!(ρs, ρsMat, δLs)
#     sum(energies) # return the kinetic part of the gs energy
# end

import Dagger

function densitymatrix_dagger!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks, 2)
    energies = Dagger.@mutable zeros(Float64, L)
    δLs = collect(keys(ρs))

    M = length(δLs); Ns = size(H)
    ρsMat = Dagger.@mutable zeros(ComplexF64, Ns..., M)

    function getspectrum(i_)
        k = ks[:, i_]
        spectrum(H(k); kwargs...)
    end

    # Launch spectrum calculations
    spectrum_tasks = []
    for i_ = 1:L
        push!(spectrum_tasks, (i_, Dagger.@spawn getspectrum(i_)))
    end

    # Launch addition tasks
    function addToRho!(ρsMat, energies, spectrum_k)
        (i_, (ϵs, U)) = fetch(spectrum_k)
        w = kweights[i_] 
        k = ks[:, i_] 

        fd = fermidirac.(real.(ϵs .- μ); T=T)
        phases = [fourierphase(-k, δL) for δL in δLs]

        Tullio.@tullio  M[i, j, n] := $w * fd[m] * phases[n] * U[i, m] * conj(U[j, m])
        ρsMat[:] .+= M[:]

        energies[i_] = real(w * sum(ϵs .* fermidirac(ϵs .- μ)))
        nothing
    end

    addition_taks = [Dagger.@spawn addToRho!(ρsMat, energies, spectrum_k) for spectrum_k in spectrum_tasks]
    wait.(addition_taks)

    flexibleformat!(ρs, collect(ρsMat), δLs)
    sum(collect(energies)) # return the kinetic part of the gs energy
end
# function densitymatrix_dagger!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)

#     L = size(ks, 2)
#     energies = Dagger.@mutable zeros(Float64, L)

#     for (δL, ρ0) in ρs
#         ρs[δL][:] .= 0.0
#     end

#     ρs1 = Dagger.@mutable ρs



#     results = [(Dagger.@spawn (j_, spectrum(H(ks[:, j_]); kwargs...))) for j_ = 1:N]

#     t = Dagger.spawn((bands, results) -> begin
#             for (j_, result) = results
#                 k = ks[:, j_]
#                 ϵs, U = result


#             end
#         end, bands, results)
#     wait(t)

#     @showprogress progressmin "Eigensolver... " for i_ = 1:L

#         k = ks[:, i_]
#         ϵs, U = spectrum(H(k); kwargs...) #@time
#         fd = fermidirac.(real.(ϵs .- μ); T=T)

#         Tullio.@tullio M[i, j] := fd[m] * U[i, m] * conj(U[j, m])

#         w = kweights[i_]
#         for δL in keys(ρs)
#             ρs[δL][:, :] .+= w .* (M .* fourierphase(-k, δL))
#         end

#         energies[i_] = real(w * sum(ϵs .* fermidirac(real(ϵs .- μ))))
#     end

#     sum(energies) # return the kinetic part of the gs energy
# end


function densitymatrix_serial!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)

    L = size(ks, 2)
    energies = zeros(Float64, L)

    for (δL, ρ0) in ρs
        ρs[δL][:] .= 0.0
    end

    @showprogress progressmin "Eigensolver... " for i_ = 1:L
        k = ks[:, i_]
        ϵs, U = spectrum(H(k); kwargs...) #@time
        fd = fermidirac.(real.(ϵs .- μ); T=T)

        Tullio.@tullio M[i, j] := fd[m] * U[i, m] * conj(U[j, m])
        # TensorOperations.@tensoropt M[i, j] := fd[m] * U[i, m] * conj(U[j, m])

        w = kweights[i_]
        for δL in keys(ρs)
            ρs[δL][:, :] .+= w .* (M .* fourierphase(-k, δL))
        end

        energies[i_] = real(w * sum(ϵs .* fermidirac(real(ϵs .- μ))))
    end

    sum(energies) # return the kinetic part of the gs energy
end


################################################################################
# KEPT FOR FUTURE REFERENCE:
################################################################################

### this is tested and working but does not really perform well:
# function densitymatrix_multithread!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64 = 0.0; T::Real = 0.01, progressmin::Int = 20, kwargs...)
#     L = size(ks, 2)

#     energies = zeros(Float64, L)
#     function spectrumf(k)
#         spectrum(H(k); kwargs...)
#     end

#     for (δL, ρ0) in ρs
#         ρs[δL][:] .= 0 #convert(SharedArray, zero(ρ0))[:]
#     end

#     Msize = size(first(values(ρs)))
#     Mtype = eltype(first(values(ρs)))

#     lk = Threads.ReentrantLock()
#     Threads.@threads for i_ in ProgressBar(1:L)#i_=1:L
#         k = ks[:, i_]
#         ϵs, U = spectrumf(k) #@time

#         M = zeros(Mtype, Msize)
#         densitymatrix!(M, ϵs .- μ, U; T = T)

#         lock(lk) do
#             for δL in keys(ρs)
#                 ρs[δL][:, :] .+= (M .* fourierphase(-k, δL))
#             end
#         end

#         energies[i_] = groundstate_sumk(real(ϵs), μ)
#     end

#     for δL in keys(ρs)
#         ρs[δL][:] ./= L
#     end

#     sum(energies) / L # return the groundstate energy
# end

# function getdensitymatrix(H, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
#     d = dim(H, ks)
#     ρ0 = zeros(ComplexF64, d, d)
#     densitymatrix!(ρ0, H, ks, μ; kwargs...)

#     ρ0
# end

# function densitymatrix!(ρ0::AbstractMatrix, δL::AbstractVector, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; T=0.01, kwargs...)

#     phase = fourierphase(-k, δL)
#     densitymatrix!(ρ0, ϵs, U; φk=phase, T=T, kwargs...)
# end

# function mapdensitymatrix!(ρ0::AbstractMatrix, δLs::AbstractVector{AbstractVector}, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; T=0.01, kwargs...)

#     M = zero(ρ0)
#     densitymatrix!(M, ϵs, U; φk=phase, T=T, kwargs...)

#     for δL=δLs 
#         ρ0[:,:] .+= M .* fourierphase(-k, δL)
#     end

#     ρ0
# end

# function densitymatrix!(ρs::AbstractHops, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; kwargs...)

#     for δL=keys(ρs)
#         densitymatrix!(ρs[δL], δL, k, ϵs, U; kwargs...)
#     end
#     ρs
# end

#### Nice idea but less performant... large overhead for copying data between processes
# function densitymatrix_pmap_async!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
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

# function densitymatrix_pmap!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
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
#         ρs[L] .= ρsnew[L]./L
#     end

#     energies # return the groundstate energy
# end

# function densitymatrix_pmap!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
#     L = size(ks,2)

#     function spectrumf(k)
#         spectrum(H(k); kwargs...)
#     end

#     zeromat, δLs = efficientzero(ρs)

#     ρ0 = convert(SharedArray, zeromat)

#     energy = sum(pmap(1:L) do i_
#         k = ks[:,i_]
#         ϵs, U = spectrumf(k) #@time

#         M = zero(ρ0[:,:,1])

#         densitymatrix!(M, ϵs.-μ, U; T=T)

#         for (j_,δL)=enumerate(δLs)
#             ρ0[:,:,j_] .+= (M .* fourierphase(-k, δL))
#         end

#         groundstate_sumk(real(ϵs), μ)
#     end)/L

#     flexibleformat!(ρs, ρ0/L, δLs)

#     energy # return the groundstate energy
# end