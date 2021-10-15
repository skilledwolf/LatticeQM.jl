
using LinearAlgebra: Hermitian

energies(H, args...; kwargs...) = geteigvals(H, args...; kwargs...)
wavefunctions(H, args...; kwargs...) = geteigvecs(H, args...; kwargs...)
spectrum(H, args...; kwargs...) = geteigen(H,args...; kwargs...)

using Distributed
using ProgressMeter
# using ProgressBars

dim(A::AbstractMatrix, x=0) = size(A,1)
dim(f::Function, x::Number) = size(f(x), 1)
dim(f::Function, x::AbstractVector) = size(f(first(x)), 1)
dim(f::Function, x::AbstractMatrix) = size(f(first(eachcol(x))), 1)

handleprojector(projector) = isa(projector, AbstractVector) ? [expvalf(p) for p in projector] : [expvalf(projector)]

"""
    bandmatrix(H, ks::Matrix{Float} [, As]; kwargs...)

Calculates the energies for operator `H(k)` for each column vector `k` of matrix `ks`.
If operators `As=[A1, A2, ...]` are given, their expectaction values are
calculated and stored for each eigenvector.

Accepts the same keywords as `geteigvals`, `geteigvecs`, `geteigen`.
In particular: `format` (`:sparse` or `:dense`) and `num_bands::Int`.

Returns a matrix where each column contains the energies for each column `k` in `ks`.
If `As` were given, a second matrix of the same format is returned containing expectation values.

### Example
```julia
using LatticeQM

lat = Geometries.honeycomb()
h = Operators.graphene(lat)
ks = kpath(lat; num_points=200)
valley = Operators.valleyoperator(lat)

bands, obs = bandmatrix(h, ks.points, valley)

```
"""
function bandmatrix(args...; multimode=:distributed, kwargs...)
    if multimode == :distributed && nprocs()>1
        bandmatrix_distributed(args...; kwargs...) # or bandmatrix_pmap ?
    # elseif multimode == :distributed2 && nprocs()>1
    #     bandmatrix_distributed(args...; kwargs...)
    elseif multimode == :multithread && Threads.nthreads()>1
        bandmatrix_multithread(args...; kwargs...)
    else
        bandmatrix_serial(args...; kwargs...)
    end
end


function bandmatrix_serial(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

    N = size(ks,2) # no. of k points
    bands = zeros(Float64, D, N)

    function energiesf(k)
        energies(H(k); kwargs...)
    end

    @showprogress (hidebar ? 10^6 : 20) "Computing bands... " for j_=1:N
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
    end

    bands
end

function bandmatrix_distributed(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

    N = size(ks,2) # no. of k points
    bands = convert(SharedArray, zeros(Float64, D, N))

    function energiesf(k)
        energies(H(k); kwargs...)
    end

    @sync @showprogress (hidebar ? 10^6 : 20) "Computing bands... "  @distributed for j_=1:N
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
    end

    convert(Array, bands)
end

function bandmatrix_pmap(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)

    bands = hcat((@showprogress (hidebar ? 10^6 : 20) "Computing bands... " pmap(x->real(energies(H(x); kwargs...)), eachcol(ks)))...)

    bands
end


function bandmatrix_multithread(H, ks; num_bands::Int=0, kwargs...)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

    N = size(ks,2) # no. of k points
    bands = zeros(Float64, D, N)

    function energiesf(k)
        energies(H(k); kwargs...)
    end

    # lk = Threads.ReentrantLock()
    Threads.@threads for j_=1:N
        # en = real.(energiesf(ks[:,j_]))
        # lock(lk) do 
        bands[:,j_] .= real.(energiesf(ks[:,j_]))
        # end
    end

    bands
end

function bandmatrix_serial(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = zeros(Float64, D, N)
    obs   = zeros(Float64, D, N, L)

    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    @showprogress (hidebar ? 10^6 : 20) "Computing bands... " for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    bands, obs
end

function bandmatrix_distributed(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = convert(SharedArray, Matrix{Float64}(undef, (D,N)))
    obs   = convert(SharedArray, Array{Float64}(undef, (D,N,L)))

    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    @sync @showprogress (hidebar ? 10^6 : 20) "Computing bands... " @distributed for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    Array(bands), Array(obs)
end

function bandmatrix_pmap(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)

    res = (@showprogress (hidebar ? 10^6 : 20) "Computing bands... " pmap(eachcol(ks)) do k
        ϵs, U = spectrum(H(k); kwargs...)
        colors = [ P(k, psi, e) for (e,psi) in zip(ϵs, eachcol(U)), P in projector ]

        ϵs, colors
    end)

    bands = hcat((x[1] for x=res)...)
    obs = cat((x[2] for x=res)...; dims=3)

    real.(bands), permutedims(obs, [1,3,2])
end

# function bandmatrix_pmap(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
#     projector = handleprojector(projector)
#     kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
#     D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

#     N = size(ks, 2) # number of k points
#     L = length(projector)
#     bands = convert(SharedArray, Matrix{Float64}(undef, (D,N)))
#     obs   = convert(SharedArray, Array{Float64}(undef, (D,N,L)))

#     function spectrumf(k)
#         spectrum(H(k); kwargs...)
#     end

#     @showprogress (hidebar ? 10^6 : 20) "Computing bands... " pmap(1:N) do j_
#         ϵs, U = spectrumf(ks[:,j_])
#         bands[:,j_] .= real.(ϵs)

#         for i_=1:size(U,2), n_=1:L
#             obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
#         end
#         nothing
#     end

#     Array(bands), Array(obs)
# end


function bandmatrix_multithread(H, ks, projector; num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = zeros(D, N)
    obs   = zeros(D, N, L)

    function spectrumf(k)
        spectrum(H(k); kwargs...)
    end

    # lk = Threads.ReentrantLock()
    Threads.@threads for j_=1:N
        ϵs, U = spectrumf(ks[:,j_])
        
        # lock(lk) do 
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
        # end
    end

    bands, obs
end


import ..Structure.Paths: DiscretePath

"""
    getbands(H, ks::DiscretePath [, As]; kwargs...)

Calculates the bands for operator `H` along discrete path `ks` and
if operators `As=[A1, A2, ...]` are given, their expectaction values are
calculated and stored for each eigenvector.

Note that ks is a discrete path object as returned by `kpath(lat::Lattice,...)`.

Accepts the same keywords as `geteigvals`, `geteigvecs`, `geteigen`.
In particular: `format` (`:sparse` or `:dense`) and `num_bands::Int`.

Returns a `BandData` object (with fields `bands`, `obs`, `path`).

### Example
```julia
using LatticeQM

lat = Geometries.honeycomb()
h = Operators.graphene(lat)
ks = kpath(lat; num_points=200)
valley = Operators.valleyoperator(lat)

bands = getbands(h, ks, valley)

using Plots
plot(bands)
```
"""
function getbands(H, ks::DiscretePath; kwargs...)
    bands = bandmatrix(H, ks; kwargs...)
    obs = nothing
    BandData(bands, obs, ks)
end

function getbands(H, ks::DiscretePath, projector; kwargs...)
    bands, obs = bandmatrix(H, ks, projector; kwargs...)
    BandData(bands, obs, ks)
end