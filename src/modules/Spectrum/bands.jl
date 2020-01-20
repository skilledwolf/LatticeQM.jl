using Distributed
using ProgressMeter

using ..Structure.Paths: DiscretePath, points
using ..TightBinding: expvalf, AnyHops, dim

handleprojector(projector::Nothing) = nothing
handleprojector(projector::Function) = [projector]
handleprojector(projector::AbstractMatrix) = handleprojector([projector])
handleprojector(projector::AnyHops) = handleprojector([projector])
function handleprojector(projector::AbstractVector)
    if length(projector) == 0
        return nothing
    end

    handle(p::Function) = p
    handle(p::AbstractMatrix) =  expvalf(p)
    handle(p::AnyHops) =  expvalf(p)

    return [handle(p) for p in projector]
end

function bandmatrix(H, ks; num_bands::Int=0, kwargs...)
    ks = points(ks)
    if !(num_bands>0)
        num_bands = dim(H, ks)
    end
    N = size(ks,2) # no. of k points
    bands = convert(SharedArray, zeros(Float64, num_bands, N))

    energiesf = energies(H; num_bands=num_bands, kwargs...)

    @sync @showprogress 1 "Computing bands..."  @distributed for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
    end

    convert(Array, bands)
end

function bandmatrix(H, ks, projector; num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    ks = points(ks)
    if !(num_bands>0)
        num_bands = dim(H, ks)
    end

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = convert(SharedArray, zeros(Float64, num_bands, N))
    obs   = convert(SharedArray, zeros(Float64, num_bands, N, L))

    spectrumf = spectrum(H; num_bands=num_bands, kwargs...)

    @sync @showprogress 1 "Computing bands..." @distributed for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    Array(bands), Array(obs)
end



function getbands(H, ks::DiscretePath; kwargs...)
    bands = bandmatrix(H, points(ks); kwargs...)
    obs = nothing
    BandData(bands, obs, ks)
end

function getbands(H, ks::DiscretePath, projector; kwargs...)
    bands, obs = bandmatrix(H, points(ks), projector; kwargs...)
    BandData(bands, obs, ks)
end

export get_bands
@legacyalias getbands get_bands


## The following code snippet might come in handy some day:
# p = Progress(nk)
# channel = RemoteChannel(()->Channel(nk), 1)
#
# @sync begin
#     # this task prints the progress bar
#     @async begin
#         io = open("bandstest.out", "w")
#         done = 0
#         while done < nk
#             (j, energies) = take!(channel)
#             writedlm(io, [j..., energies...]')
#             bands[:,j] = energies
#
#             next!(p)
#             done = done + 1
#         end
#         close(io)
#     end
#
#
#     # this task does the computation
#     @async begin
#         @sync @distributed for j_=1:nk
#             k=kgrid[:,j_]
#
#             ϵs, U = Σ(k)
#
#             put!(channel, (j_, real.(ϵs)))
#
#         end
#     end
# end