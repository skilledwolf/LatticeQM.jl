using Distributed
using ProgressMeter

const energies = geteigvals
const wavefunctions = geteigvecs
const spectrum = geteigen

const IMAG_THRESHOLD = 1e-7::Float64
const PROGRESSBAR_MINTIME = 1::Int
const PROGRESSBAR_DIAG_DEFAULTLABEL = "Diagonalization"::String
const PROGRESSBAR_SHOWDEFAULT = true::Bool
const PROGRESSBAR_BANDMATRIX_DEFAULTKWARGS = Dict(:hidebar => !PROGRESSBAR_SHOWDEFAULT, :progress_label => PROGRESSBAR_DIAG_DEFAULTLABEL)

################################################################################
# dim helper functions (Should maybe be moved to Utils?)
################################################################################

dim(A::AbstractMatrix, x=0) = size(A,1)
dim(f::Function, x::Number) = size(f(x), 1)
dim(f::Function, x::AbstractVector) = size(f(first(x)), 1)
dim(f::Function, x::AbstractMatrix) = size(f(first(eachcol(x))), 1)

################################################################################
# Expectation value helper type/functions
################################################################################

abstract type AbstractExpvalMap end
struct ExpvalMap{T} <: AbstractExpvalMap
    op::T
end
(O::ExpvalMap{<:AbstractMatrix})(_, ψ, _) = real.(dot(ψ, O.op, ψ))
(O::ExpvalMap)(k, ψ, _) = real.(dot(ψ, O.op(k), ψ))
ExpvalMap(op::Function) = ExpvalMap{Function}(op)
ExpvalMap(op::T) where T<:AbstractMatrix = ExpvalMap{T}(op)
ExpvalMap(op::ExpvalMap) = op

handleprojector() = AbstractExpvalMap[]
handleprojector(projector::AbstractExpvalMap) = typeof(projector)[projector]
handleprojector(projector) = handleprojector([projector]) # Anything that is not a vector is converted to a vector first
handleprojector(projector::AbstractVector) = [ExpvalMap(op) for op in projector]
handleprojector(projector::AbstractVector{<:AbstractExpvalMap}) = projector

################################################################################
# Low-level helper functions
################################################################################

function assert_realeigvals(ϵs)
    imag_check = imag.(ϵs) .< IMAG_THRESHOLD
    @assert all(imag_check) "Imaginary eigenvalues encountered!: $(ϵs[.!imag_check])"
end

function bandmatrix_size(H, ks; kwargs...)
    D = get(kwargs, :num_bands, dim(H, ks))::Int # check if num_bands is given
    N = size(ks, 2)
    D, N
end

function compute_bandexpvals!(obs::AbstractMatrix, projectors::AbstractVector, k, ϵs, U)
    for i_ = axes(U, 2), n_ = eachindex(projectors)
        obs[i_, n_] = projectors[n_](k, U[:, i_], ϵs[i_])
    end
    obs
end

insertbands!(bands::AbstractVector, ϵs::AbstractVector) = (assert_realeigvals(ϵs); bands .= real.(ϵs); bands)

function insertbands_bandexpvals_k!(bands::AbstractVector, obs::AbstractMatrix, H::AbstractMatrix, k, projectors; kwargs...)
    ϵs, U = spectrum(H; kwargs...)
    assert_realeigvals(ϵs)
    bands .= real.(ϵs)
    compute_bandexpvals!(obs, projectors, k, ϵs, U)
    bands, obs
end

function bandmatrix_preallocate(H, ks, projectors; kwargs...)
    D, N = bandmatrix_size(H, ks; kwargs...)
    L = length(projectors)
    bands = zeros(Float64, D, N)
    obs = zeros(Float64, D, N, L)
    bands, obs
end

################################################################################
# Main functions
################################################################################

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
function bandmatrix(H, ks, projectors...; multimode=:distributed, kwargs...)
    ks = Structure.points(ks) # sanatize ks to be a matrix
    projectors = handleprojector(projectors...) # sanatize projectors
    bands, obs = bandmatrix_preallocate(H, ks, projectors; kwargs...)

    if multimode == :distributed && nprocs() > 1
        bands, obs = SharedArray(bands), SharedArray(obs) # convert to shared arrays
        return bandmatrix_distributed!(bands, obs, H, ks, projectors; kwargs...)
    if multimode == :multithreaded && Threads.nthreads() > 1
        return bandmatrix_multithreaded!(bands, obs, H, ks, projectors; kwargs...)
    else 
        return bandmatrix_serial!(bands, obs, H, ks, projectors; kwargs...)
    end
end

function bandmatrix_serial!(bands, obs, H, ks, projectors; hidebar=!PROGRESSBAR_SHOWDEFAULT, progress_label=PROGRESSBAR_DIAG_DEFAULTLABEL, kwargs...)
    @showprogress dt = PROGRESSBAR_MINTIME desc = progress_label enabled = !hidebar for j_ = axis(bands, 2)
        @views insertbands_bandexpvals_k!(bands[:, j_], obs[:, j_, :], H(ks[:, j_]), ks[:, j_], projectors; kwargs...)
    end
    bands, obs
end

function bandmatrix_distributed!(bands, obs, H, ks, projectors; hidebar=!PROGRESSBAR_SHOWDEFAULT, progress_label=PROGRESSBAR_DIAG_DEFAULTLABEL, kwargs...)
    @sync @showprogress dt=PROGRESSBAR_MINTIME desc=progress_label enabled=!hidebar @distributed for j_ = axis(bands, 2)
        @views insertbands_bandexpvals_k!(bands[:, j_], obs[:, j_, :], H(ks[:, j_]), ks[:, j_], projectors; kwargs...)
    end
    bands, obs
end

function bandmatrix_multithreaded!(bands, obs, H, ks, projectors; hidebar=!PROGRESSBAR_SHOWDEFAULT, progress_label=PROGRESSBAR_DIAG_DEFAULTLABEL, kwargs...)
    p = Progress(size(bands, 2), PROGRESSBAR_MINTIME, progress_label; enabled = !hidebar)
    @showprogress Threads.@threads for j_ = axes(bands, 2)
        @views insertbands_bandexpvals_k!(bands[:, j_], obs[:, j_, :], H(ks[:, j_]), ks[:, j_], projectors; kwargs...)
        next!(p)
    end
    finish!(p)
    bands, obs
end

################################################################################
# Main public user interface 
################################################################################

# import ..Structure.Paths: DiscretePath

"""
    getbands(H, ks::Union{DiscretePath, AbstractMatrix} [, As]; kwargs...)

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
function getbands(H, ks, projectors...; kwargs...)
    bands, obs = bandmatrix(H, ks, projectors...; kwargs...)
    BandData(bands, obs, ks)
end


################################################################################
# Dagger parallelization
################################################################################

# import Dagger
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