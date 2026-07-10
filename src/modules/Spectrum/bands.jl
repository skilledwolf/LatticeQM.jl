using Distributed
using SharedArrays
using ProgressMeter
import SparseArrays
import LatticeQM.Eigen
import LatticeQM.Parallel

# Build H(k) into a preallocated buffer. Generic fallback: materialise `H(k)`
# (which may allocate inside the user-supplied operator) and copy. The
# `AbstractHops` specialisation in TightBinding/types.jl uses `fouriersum!`
# for a true zero-allocation path — that's what removes ~1 GB allocation
# churn / ~30% GC overhead on a TBG-N=5 bandmatrix call.
@inline _build_H!(out::AbstractMatrix, H, k) = (out .= H(k); out)

# Bloch-matrix scratch utilities, used by every k-loop (bandmatrix, dos,
# density, density-matrix, Green, optical conductivity).
#
# `bloch_buffer(H, ks; format)` returns a buffer for in-place Hk construction:
#   - `format=:dense`  → preallocated dense Matrix{ComplexF64}, reused per k
#   - `format=:sparse` → `nothing` (sparse H(k) is built fresh per k; the
#     per-k alloc is O(nnz), much smaller than O(N²), and crucially keeps
#     the matrix sparse so the eigensolver does sparse LU during shift-invert)
#
# `bloch!(buf, H, k)` returns the Bloch matrix at k. For a dense buffer this
# is zero-allocation when H is an `AbstractHops` (uses `fouriersum!`). For
# `nothing` it just calls `H(k)`.
bloch_buffer(H, ks; format=:dense, kwargs...) =
    format === :sparse ? nothing : bandmatrix_Hcache(H, ks)

@inline bloch!(buf::Nothing, H, k) = H(k)
@inline bloch!(buf::AbstractMatrix, H, k) = (_build_H!(buf, H, k); buf)

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
    imag_check = abs.(imag.(ϵs)) .< IMAG_THRESHOLD   # one-sided check let large negative imag parts through
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
If operators `As=[A1, A2, ...]` are given, their expectation values are
calculated and stored for each eigenvector.

Accepts the same keywords as `geteigvals`, `geteigvecs`, `geteigen`. In particular:
`format` (`:sparse` or `:dense`) and `num_bands::Int`.

Parallelism is selected via `multimode`:
- `:serial` — no parallelism
- `:multithreaded` / `:threaded` — `Threads.@spawn` over the default pool
- `:distributed` — `pmap` over `Distributed.workers()`
- `:auto` (default) — distributed if workers exist, else threads if >1, else serial

For full control, pass a `Parallel.Executor` directly via the `executor` kwarg
(e.g. `executor=Parallel.ThreadedExec(8; schedule=:static)`).

Returns `(bands, obs)`: bands is `num_bands × num_k`, obs is the matching tensor
of projector expectation values.

### Example
```julia
using LatticeQM
lat = Geometries.honeycomb()
h = Operators.graphene(lat)
ks = kpath(lat; num_points=200)
bands, obs = bandmatrix(h, ks.points, Operators.valley(lat))
```
"""
function bandmatrix(H, ks, projectors...;
                    multimode=:auto,
                    executor::Union{Nothing,Parallel.Executor}=nothing,
                    progress_label=PROGRESSBAR_DIAG_DEFAULTLABEL,
                    hidebar=!PROGRESSBAR_SHOWDEFAULT,
                    kwargs...)

    ks = Structure.points(ks)
    projectors = handleprojector(projectors...)
    bands, obs = bandmatrix_preallocate(H, ks, projectors; kwargs...)

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor

    # Distributed shared-memory optimisation: convert H to SharedDenseHops
    # so workers don't each receive a full copy via serialisation.
    if exec isa Parallel.DistributedExec
        H = sanatize_distributed_hamiltonian(H)
        bands = SharedArray(bands)
        obs = SharedArray(obs)
    end

    Parallel.configure_blas!(exec; verbose=false)

    progressbar = Progress(size(ks, 2);
                           dt=PROGRESSBAR_MINTIME,
                           desc=progress_label,
                           enabled=!hidebar)

    Parallel.kspace_foreach!(ks, exec;
        scratch_factory = () -> bandmatrix_scratch(H, ks; kwargs...),
        progress = progressbar) do scratch, j, k
        _band_kpoint!(scratch, j, k, bands, obs, H, projectors; kwargs...)
    end
    finish!(progressbar)

    if exec isa Parallel.DistributedExec
        bands, obs = sdata(bands), sdata(obs)
    end

    bands, obs
end

# Per-task scratch: H buffer reused across every k handled by the same task.
# Delegates to `bloch_buffer` / `bloch!` so dos / density / green / etc. can
# all use the same per-task pattern.
bandmatrix_scratch(H, ks; format=:dense, kwargs...) =
    (Hcache = bloch_buffer(H, ks; format=format),)

function _band_kpoint!(scratch, j, k, bands, obs, H, projectors; kwargs...)
    Hk = bloch!(scratch.Hcache, H, k)
    if isempty(projectors)
        # Eigenvectors aren't needed (no observables / projectors). Skipping
        # them via `geteigvals!` instead of `geteigen!` is ~2.5× faster on
        # the dense Hermitian path (LAPACK `heevr` only does the
        # eigenvalue half). The chempot solve hits this branch on every
        # call.
        energies_k = Eigen.geteigvals!(Hk; kwargs...)
        @views insertbands!(bands[:, j], energies_k)
    else
        energies_k, U = Eigen.geteigen!(Hk; kwargs...)
        @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], energies_k, U, k, projectors)
    end
    return nothing
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
