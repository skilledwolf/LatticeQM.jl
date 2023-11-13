
using LinearAlgebra: Hermitian

const energies = geteigvals
const wavefunctions = geteigvecs
const spectrum = geteigen

using Distributed
using ProgressMeter
# using ProgressBars

dim(A::AbstractMatrix, x=0) = size(A,1)
dim(f::Function, x::Number) = size(f(x), 1)
dim(f::Function, x::AbstractVector) = size(f(first(x)), 1)
dim(f::Function, x::AbstractMatrix) = size(f(first(eachcol(x))), 1)

function expvalf(𝑶::AbstractMatrix)
    f(k, ψ, ϵ) = real.(dot(ψ, 𝑶, ψ))
    f
end

function expvalf(𝑶::Function)
    f(k, ψ, ϵ) = real.(dot(ψ, 𝑶(k), ψ))
    f
end

handleprojector(projector) = isa(projector, AbstractVector) ? map(expvalf, projector) : [expvalf(projector)]

import ..Structure

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
function bandmatrix(H, ks, args...; multimode = :distributed, kwargs...)
    ks = Structure.points(ks)

    return bandmatrix_distributed(H, ks, args...; kwargs...)

    if multimode == :distributed && nprocs() > 1
        bandmatrix_distributed(H, ks, args...; kwargs...) # or bandmatrix_pmap ?
    # elseif multimode == :distributed2 && nprocs()>1
    #     bandmatrix_distributed(H, ks, args...; kwargs...)
    elseif multimode == :multithread && Threads.nthreads() > 1
        bandmatrix_multithread(H, ks, args...; kwargs...)
    else
        bandmatrix_serial(H, ks, args...; kwargs...)
    end
end


function bandmatrix_serial(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands>0) ? num_bands : dim(H,ks) # matrix dimension

    N = size(ks,2) # no. of k points
    bands = zeros(Float64, D, N)

    @showprogress (hidebar ? 10^6 : 20) "Computing bands... " for j_=1:N
        bands[:,j_] .= energies(H(ks[:,j_]); kwargs...)
    end

    bands
end

import Dagger
# function bandmatrix_distributed(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
#     kwargs = num_bands == 0 ? kwargs : Dict(kwargs..., :num_bands => num_bands)
#     D = (num_bands > 0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

#     N = size(ks, 2)
#     bands = Dagger.@mutable zeros(Float64, D, N)

#     results = [(Dagger.@spawn (j_, energies(H(ks[:, j_]); kwargs...))) for j_ = 1:N]
#     t = Dagger.spawn((bands, results) -> begin
#             for (j_, result) = results
#                 bands[:, j_] .= real.(fetch(result))
#             end
#         end, bands, results)
#     wait(t)

#     collect(bands)
# end
# function bandmatrix_distributed(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
#     kwargs = num_bands == 0 ? kwargs : Dict(kwargs..., :num_bands => num_bands)
#     D = (num_bands > 0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

#     N = size(ks, 2) # no. of k points
#     bands = Dagger.@mutable zeros(Float64, D, N)

#     # @time energies(H(ks[:, 2]))

#     # @time bands = begin
#     #     ts = []
#     #     for j_ = 1:N
#     #         hk = Dagger.@spawn H(ks[:, j_])
#     #         en = Dagger.@spawn energies(hk; kwargs...)
            
#     #         t = Dagger.spawn((bands, energies, j_) -> begin
#     #             bands[:, j_] .= real.(energies)
#     #             nothing
#     #         end, bands, en, j_)
#     #         push!(ts, t)
#     #     end
#     #     wait.(ts)
#     #     collect(bands)
#     # end

#     # @time begin
#     #     results = [(Dagger.@spawn (j_, energies(H(ks[:, j_]); kwargs...))) for j_ = 1:N]
#     #     # ts = [Dagger.spawn((bands, energies, j_) -> begin
#     #     #     bands[:, j_] .= real.(energies)
#     #     #     nothing
#     #     # end, bands, result, j_) for (j_, result)=results]
#     #     # wait.(ts)
#     #     t = Dagger.spawn((bands, results) -> begin
#     #             for (j_, result) = results
#     #                 bands[:, j_] .= real.(fetch(result))
#     #             end
#     #         end, bands, results)
#     #     wait(t)
#     #     bands = collect(bands)
#     # end


#     results = [(Dagger.@spawn (j_, energies(H(ks[:,j_]); kwargs...))) for j_=1:N]
#     t = Dagger.spawn((bands, results) -> begin
#         for (j_,result) = results 
#             bands[:, j_] .= real.(fetch(result))
#         end
#     end, bands, results)
#     wait(t)
#     bands = collect(bands)

#     # println("norm ", LinearAlgebra.norm(bands))
#     # println("BANDMATRIX DONE")
#     # println("Stopping workers ...")
#     # println(workers())
#     # rmprocs(workers())
#     # exit()

#     bands
# end

function bandmatrix_distributed(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands > 0) ? num_bands : dim(H, ks) # matrix dimension

    N = size(ks,2) # no. of k points
    bands = SharedArray(Array{Float64}(undef, D, N))

    @sync @showprogress (hidebar ? 10^6 : 5) "Computing bands... "  @distributed for j_=1:N
        bands[:, j_] .= real.(energies(H(ks[:, j_]); kwargs...))
    end

    # sdata(bands)
    bands
end

# function bandmatrix_pmap(H, ks; hidebar=false, num_bands::Int=0, kwargs...)
#     kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)

#     bands = hcat((@showprogress (hidebar ? 10^6 : 20) "Computing bands... " pmap(x->real(energies(H(x); kwargs...)), eachcol(ks)))...)

#     bands
# end


function bandmatrix_serial(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands > 0) ? num_bands : dim(H, ks) # matrix dimension

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = zeros(Float64, D, N)
    obs   = zeros(Float64, D, N, L)

    function spectrumf(k)
        
    end

    @showprogress (hidebar ? 10^6 : 20) "Computing bands... " for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrum(H(ks[:, j_]); kwargs...)
        @assert all(imag.(ϵs) .< 1e-10) "Imaginary eigenvalues encountered!"
        bands[:,j_] .= real.(ϵs)

        for i_=1:size(U,2), n_=1:L
            obs[i_,j_,n_] = projector[n_](ks[:,j_],U[:,i_],ϵs[i_])
        end
    end

    bands, obs
end

# function bandmatrix_distributed(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
#     projector = handleprojector(projector)
#     kwargs = num_bands == 0 ? kwargs : Dict(kwargs..., :num_bands => num_bands)
#     D = (num_bands > 0) ? num_bands : size(H(first(eachcol(ks))), 1) # matrix dimension

#     N = size(ks, 2) # no. of k points
#     L = length(projector)
#     bands = Dagger.@mutable zeros(Float64, D, N)
#     obs = Dagger.@mutable zeros(Float64, D, N, L)

#     results = [(Dagger.@spawn (j_, spectrum(H(ks[:, j_]); kwargs...))) for j_ = 1:N]
#     t = Dagger.spawn((bands, obs, results) -> begin
#             for (j_, result) = results
#                 ϵs, U = fetch(result)
#                 bands[:, j_] .= real.(ϵs)

#                 for i_ = 1:size(U, 2), n_ = 1:L
#                     obs[i_, j_, n_] = projector[n_](ks[:, j_], U[:, i_], ϵs[i_])
#                 end
#             end
#         end, bands, obs, results)
#     wait(t)

#     collect(bands), collect(obs)
# end

function bandmatrix_distributed(H, ks, projector; hidebar=false, num_bands::Int=0, kwargs...)
    projector = handleprojector(projector)
    kwargs = num_bands==0 ? kwargs : Dict(kwargs..., :num_bands=>num_bands)
    D = (num_bands > 0) ? num_bands : dim(H, ks) # matrix dimension

    N = size(ks, 2) # number of k points
    L = length(projector)
    bands = convert(SharedArray, Matrix{Float64}(undef, (D,N)))
    obs   = convert(SharedArray, Array{Float64}(undef, (D,N,L)))

    @sync @showprogress (hidebar ? 10^6 : 20) "Computing bands... " @distributed for j_=1:N
#     @showprogress 1 "Computing bands..." for j_=1:N
        ϵs, U = spectrum(H(ks[:,j_]); kwargs...)
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