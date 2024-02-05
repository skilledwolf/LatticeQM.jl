using Distributed
using SharedArrays
using ProgressMeter
import SparseArrays
import LatticeQM.Eigen

const IMAG_THRESHOLD = 1e-9::Float64
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

function compute_chunks(total_length, num_chunks)
    chunk_size = ceil(Int, total_length / num_chunks)
    return [(1+(i-1)*chunk_size):min(i * chunk_size, total_length) for i in 1:num_chunks]
end

function bandmatrix_size(H, ks; kwargs...)
    num_bands = get(kwargs, :num_bands, dim(H, ks))::Int # check if num_bands is given
    num_k = size(ks, 2)
    num_bands, num_k
end

function bandmatrix_preallocate(H, ks, projectors; kwargs...)
    num_bands, num_k = bandmatrix_size(H, ks; kwargs...)
    num_expvals = length(projectors)
    bands = zeros(Float64, num_bands, num_k)
    obs = zeros(Float64, num_bands, num_k, num_expvals)
    bands, obs
end

function bandmatrix_Hcache(H, ks)
    D = dim(H, ks)
    zeros(ComplexF64, D, D)
end

function bandmatrix_Ucache(H, ks; kwargs...)
    D = dim(H, ks)
    num_bands, num_k = bandmatrix_size(H, ks; kwargs...)
    zeros(ComplexF64, D, num_bands, num_k)
end

function assert_realeigvals(ϵs)
    imag_check = imag.(ϵs) .< IMAG_THRESHOLD
    @assert all(imag_check) "Imaginary eigenvalues encountered!: $(ϵs[.!imag_check])"
end

function compute_bandexpvals!(obs::AbstractMatrix, projectors::AbstractVector, k, ϵs, U)
    for i_ = axes(U, 2), n_ = eachindex(projectors)
        obs[i_, n_] = projectors[n_](k, U[:, i_], ϵs[i_])
    end
    obs
end

insertbands!(bands::AbstractVector, ϵs::AbstractVector) = (assert_realeigvals(ϵs); bands .= real.(ϵs); bands)

function insertbands_bandexpvals_k!(bands::AbstractVector, obs::AbstractMatrix, ϵs, U, k, projectors)
    # ϵs, U = Eigen.geteigen(H; kwargs...)
    # ϵs, U = spectrum_k
    assert_realeigvals(ϵs)
    bands .= real.(ϵs)
    compute_bandexpvals!(obs, projectors, k, ϵs, U)
    bands, obs
end

################################################################################
# Main functions
################################################################################

import LatticeQM.Structure

sanatize_distributed_hamiltonian(H) = H

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
function bandmatrix(H, ks, projectors...; multimode=:distributed, progress_label=PROGRESSBAR_DIAG_DEFAULTLABEL, hidebar=!PROGRESSBAR_SHOWDEFAULT, kwargs...)

    ks = Structure.points(ks) # sanatize ks to be a matrix
    projectors = handleprojector(projectors...) # sanatize projectors
    bands, obs = bandmatrix_preallocate(H, ks, projectors; kwargs...)

    progressbar = Progress(size(ks, 2); dt=PROGRESSBAR_MINTIME, desc=progress_label, enabled=!hidebar)

    if multimode == :distributed && nprocs() > 1
        H = sanatize_distributed_hamiltonian(H)
        bands, obs = SharedArray(bands), SharedArray(obs) # convert to shared arrays
        # bandmatrix_distributed!(bands, obs, H, ks, projectors, progressbar; kwargs...)
        bandmatrix_pmap!(bands, obs, H, ks, projectors, progressbar; kwargs...)
        bands, obs = sdata(bands), sdata(obs) # convert back to normal arrays
    elseif multimode == :multithreaded && Threads.nthreads() > 1 && get(kwargs, :format, :dense) != :sparse #Arpack.eigs is not thread-safe
        bandmatrix_multithreaded!(bands, obs, H, ks, projectors, progressbar; kwargs...)
    else 
        bandmatrix_serial!(bands, obs, H, ks, projectors, progressbar; kwargs...)
    end
    finish!(progressbar)

    bands, obs
end

function bandmatrix_serial!(bands, obs, H, ks, projectors, progressbar=nothing; progress_channel=nothing, kwargs...)

    Hcache = bandmatrix_Hcache(H, ks) # local cache
    Ucache = bandmatrix_Ucache(H, ks) # local cache

    for j_ = axes(ks, 2)
        Hcache .= H(ks[:, j_])
        energies_k, Ucache = Eigen.geteigen!(Hcache; kwargs...) # Hermitian(Hcache)
        @views insertbands_bandexpvals_k!(bands[:, j_], obs[:, j_, :], energies_k, Ucache, ks[:, j_], projectors)

        !isnothing(progress_channel) && put!(progress_channel, true)
        !isnothing(progressbar) && ProgressMeter.next!(progressbar)
    end
    Hcache = nothing
    Ucache = nothing
    GC.gc() # eagerly collect garbage, especially for the local caches
    bands, obs
end

function bandmatrix_distributed!(bands::SharedArray, obs::SharedArray, H, ks, projectors, progressbar::ProgressMeter.AbstractProgress; kwargs...)
    channel = RemoteChannel(() -> Channel{Bool}(), 1)
    @sync begin
        @async while take!(channel)
            ProgressMeter.next!(progressbar)
        end

        @async begin
            bandmatrix_distributed!(bands, obs, H, ks, projectors; progress_channel=channel, kwargs...)
            put!(channel, false)
        end
    end
    bands, obs
end

function bandmatrix_distributed!(bands::SharedArray, obs::SharedArray, H, ks, projectors; num_chunks=:auto, kwargs...)
    nks = size(ks, 2)
    num_chunks = num_chunks == :auto ? nworkers() : num_chunks
    chunks = compute_chunks(nks, num_chunks)

    @sync @distributed for j_=1:num_chunks
        chunk = chunks[j_]
        @views bandmatrix_serial!(bands[:, chunk], obs[:, chunk, :], H, ks[:, chunk], projectors; kwargs...)
    end
    GC.gc()

    bands, obs
end

function bandmatrix_pmap!(bands::SharedArray, obs::SharedArray, H, ks, projectors, progressbar; num_chunks=:auto, kwargs...)
    nks = size(ks, 2)
    num_chunks = num_chunks == :auto ? nworkers() : num_chunks
    chunks = compute_chunks(nks, num_chunks)

    progress_pmap(chunks, progress=progressbar) do chunk
        @views bandmatrix_serial!(bands[:, chunk], obs[:, chunk, :], H, ks[:, chunk], projectors; kwargs...)
    end
    GC.gc()

    bands, obs
end

function bandmatrix_multithreaded!(bands::AbstractArray, obs::AbstractArray, H, ks, projectors, progressbar; num_chunks=:auto, kwargs...)
    nks = size(ks, 2)
    num_chunks = num_chunks == :auto ? Threads.nthreads() : num_chunks
    chunks = compute_chunks(nks, num_chunks)

    Threads.@threads for j_ = 1:num_chunks
        chunk = chunks[j_]
        @views bandmatrix_serial!(bands[:, chunk], obs[:, chunk, :], H, ks[:, chunk], projectors, progressbar; kwargs...)
    end
    GC.gc()

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
# Calculate bandgap
################################################################################

import LatticeQM.Structure: regulargrid

function bandgap_filling(H, filling::Real; klin=30, kwargs...)
    kgrid = regulargrid(; nk=klin^2)
    bandgap_filling(H, kgrid, filling; kwargs...)
end

function bandgap_filling(H, ks, filling::Real; multimode=:distributed, kwargs...)
    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(H, ks; multimode=multimode)[1] # dense diagonalization (default)!
    μ = chemicalpotential(bands, ks, filling; kwargs...)

    bandgap_energy(bands, μ)
end

function bandgap(H, μ::Real=0.0; klin=10, multimode=:distributed) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    kgrid = regulargrid(; nk=klin^2)
    bandgap_energy(bandmatrix(H, kgrid; multimode=multimode)[1], μ)
end

function bandgap_energy(H, ks, μ::Real=0.0, multimode=:distributed) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    bandgap_energy(bandmatrix(H, ks; multimode=multimode)[1], μ)
end

function bandgap_energy(bands::AbstractMatrix, μ::Real)
    minimum(bands[bands.>=μ]) - maximum(bands[bands.<=μ])
end
