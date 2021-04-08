using Distributed
using ProgressMeter
# using ProgressBars

import ..Structure.Paths: DiscretePath, points
import ..TightBinding: expvalf, Hops, Hops, dim

handleprojector(projector::Nothing) = nothing
handleprojector(projector::Function) = [projector]
handleprojector(projector::AbstractMatrix) = handleprojector([projector])
handleprojector(projector::Hops) = handleprojector([projector])
function handleprojector(projector::AbstractVector)
    if length(projector) == 0
        return nothing
    end

    handle(p::Function) = p
    handle(p::AbstractMatrix) =  expvalf(p)
    handle(p::Hops) =  expvalf(p)

    return [handle(p) for p in projector]
end

"""
    bandmatrix(H, ks::Matrix{Float} [, As]; kwargs...)

Calculates the energies for operator `H(k)` for each column vector `k` of matrix `ks`.
If operators `As=[A1, A2, ...]` are given, their expectaction values are
calculated and stored for each eigenvector.

Accepts the same keywords as `Algebra.geteigvals`, `Algebra.geteigvecs`, `Algebra.geteigen`.
In particular: `format` (`:sparse` or `:dense`) and `num_bands::Int`.

Returns a matrix where each column contains the energies for each column `k` in `ks`.
If `As` were given, a second matrix of the same format is returned containing expectation values.

### Example
```julia
using LatticeQM

lat = Geometries2D.honeycomb()
h = Operators.graphene(lat)
ks = kpath(lat; num_points=200)
valley = Operators.valleyoperator(lat)

bands, obs = bandmatrix(h, ks.points, valley)

```
"""
function bandmatrix(args...; multimode=:distributed, kwargs...)
    if multimode == :distributed && nprocs()>1
        bandmatrix_distributed(args...; kwargs...)
    elseif multimode == :multithread && Threads.nthreads()>1
        bandmatrix_multithread(args...; kwargs...)
    else
        bandmatrix_serial(args...; kwargs...)
    end
end

function bandmatrix_serial(H, ks; num_bands::Int=0, kwargs...)
    ks = points(ks)
    if !(num_bands>0)
        num_bands = dim(H, ks)
    end
    N = size(ks,2) # no. of k points
    bands = zeros(Float64, num_bands, N)

    energiesf = energies(H; num_bands=num_bands, kwargs...)

    @showprogress 20 "Computing bands... " for j_=1:N
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
    end

    convert(Array, bands)
end

function bandmatrix_distributed(H, ks; num_bands::Int=0, kwargs...)
    ks = points(ks)
    if !(num_bands>0)
        num_bands = dim(H, ks)
    end
    N = size(ks,2) # no. of k points
    bands = convert(SharedArray, zeros(Float64, num_bands, N))

    energiesf = energies(H; num_bands=num_bands, kwargs...)

    @sync @showprogress 20 "Computing bands... "  @distributed for j_=1:N
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
    end

    convert(Array, bands)
end


function bandmatrix_multithread(H, ks; num_bands::Int=0, kwargs...)
    ks = points(ks)
    if !(num_bands>0)
        num_bands = dim(H, ks)
    end
    N = size(ks,2) # no. of k points
    bands = zeros(Float64, num_bands, N)

    energiesf = energies(H; num_bands=num_bands, kwargs...)

    # lk = Threads.ReentrantLock()
    Threads.@threads for j_=1:N

        # en = real.(energiesf(ks[:,j_]))

        # lock(lk) do 
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
        # end
    end

    bands
end

function bandmatrix_serial(H, ks, projector; num_bands::Int=0, kwargs...)
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

    @showprogress 20 "Computing bands... " for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    Array(bands), Array(obs)
end

function bandmatrix_distributed(H, ks, projector; num_bands::Int=0, kwargs...)
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

    @sync @showprogress 20 "Computing bands... " @distributed for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    Array(bands), Array(obs)
end


function bandmatrix_multithread(H, ks, projector; num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    ks = points(ks)
    if !(num_bands>0)
        num_bands = dim(H, ks)
    end

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = zeros(num_bands, N)
    obs   = zeros(num_bands, N, L)

    spectrumf = spectrum(H; num_bands=num_bands, kwargs...)

    lk = Threads.ReentrantLock()
    Threads.@threads for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        
        lock(lk) do 
            bands[:,j_] .= real.(ϵs)

            for i_=1:size(U,2), n_=1:L
                obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
            end
        end
    end

    bands, obs
end


"""
    getbands(H, ks::DiscretePath [, As]; kwargs...)

Calculates the bands for operator `H` along discrete path `ks` and
if operators `As=[A1, A2, ...]` are given, their expectaction values are
calculated and stored for each eigenvector.

Note that ks is a discrete path object as returned by `kpath(lat::Lattice,...)`.

Accepts the same keywords as `Algebra.geteigvals`, `Algebra.geteigvecs`, `Algebra.geteigen`.
In particular: `format` (`:sparse` or `:dense`) and `num_bands::Int`.

Returns a `BandData` object (with fields `bands`, `obs`, `path`).

### Example
```julia
using LatticeQM

lat = Geometries2D.honeycomb()
h = Operators.graphene(lat)
ks = kpath(lat; num_points=200)
valley = Operators.valleyoperator(lat)

bands = getbands(h, ks, valley)

using Plots
plot(bands)
```
"""

function getbands(H, ks::DiscretePath; kwargs...)
    bands = bandmatrix(H, points(ks); kwargs...)
    obs = nothing
    BandData(bands, obs, ks)
end

function getbands(H, ks::DiscretePath, projector; kwargs...)
    bands, obs = bandmatrix(H, points(ks), projector; kwargs...)
    BandData(bands, obs, ks)
end


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