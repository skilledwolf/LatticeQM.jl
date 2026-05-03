using ProgressMeter

import LatticeQM.Utils: fermidirac
import LatticeQM.TightBinding: Hops, AbstractHops
import LatticeQM.Spectrum: dim
import LatticeQM.TightBinding: fourierphase
import LatticeQM.Structure: Mesh, meshweights
import LatticeQM.Parallel

const PROGRESSBAR_DENSMAT_DEFAULTLABEL = "Density matrix"::String

function densitymatrix!(ρ_k::AbstractMatrix, fd_k::AbstractVector, U_k::AbstractMatrix)
    mul!(ρ_k, U_k, Diagonal(fd_k) * U_k') # U_k * Diagonal(fd_k) * U_k'
    ρ_k
end

function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; kwargs...)
    L = size(ks, 2); kweights = fill(1/L, L)
    getdensitymatrix!(ρs, H, ks, kweights, μ; kwargs...)
end

function getdensitymatrix!(ρs::Hops, H, kgrid::Mesh, μ::Float64=0.0; kwargs...)
    ks = kgrid.points
    kweights = meshweights(kgrid)
    getdensitymatrix!(ρs, H, ks, kweights, μ; kwargs...)
end

function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0;
                            multimode::Symbol=:auto,
                            executor::Union{Nothing,Parallel.Executor}=nothing,
                            progress_label=PROGRESSBAR_DENSMAT_DEFAULTLABEL,
                            hidebar=!Spectrum.PROGRESSBAR_SHOWDEFAULT,
                            T::Real=0.01,
                            kwargs...)

    for δL in keys(ρs)
        ρs[δL] .= 0
    end

    progressbar = Progress(size(ks, 2);
                           dt=Spectrum.PROGRESSBAR_MINTIME,
                           desc=progress_label,
                           enabled=!hidebar)

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    out = densitymatrix_compute_add!(ρs, H, ks, kweights, μ, T, exec, progressbar; kwargs...)
    finish!(progressbar)
    out
end

###################################################################################################
###################################################################################################
###################################################################################################

import LatticeQM.Spectrum
import LatticeQM.Eigen
import LatticeQM.TightBinding
import LatticeQM.TightBinding: efficientzero, flexibleformat!, fourierphase
import LatticeQM.TightBinding: SharedDenseHops, Hops


# SharedDenseHops case: take a per-chunk view onto its underlying buffers and
# delegate to the SubarrayHops path (which is what kspace_reduce! understands).
densitymatrix_compute_add!(ρs::SharedDenseHops, args...; kwargs...) =
    densitymatrix_compute_add!(TightBinding.gethopsview(ρs), args...; kwargs...)

# Per-task scratch: a Hops accumulator with the same δL keys / matrix shapes
# as the master. Allocated once per chunk, reused across every k in the chunk.
function _densitymatrix_local(ρs::AbstractHops)
    δLs = keys(ρs)
    isempty(δLs) && return TightBinding.Hops(Dict{Vector{Int}, Matrix{ComplexF64}}())
    data = Dict{eltype(δLs), Matrix{ComplexF64}}()
    for δL in δLs
        data[δL] = zeros(ComplexF64, size(ρs[δL]))
    end
    TightBinding.Hops(data)
end

# kspace_reduce! requires `output` to support `zero(output)` and `.+=` —
# overload both for Hops so the primitive composes cleanly.
Base.zero(h::Hops) = _densitymatrix_local(h)
function _add_hops!(dst::Hops, src::Hops)
    for δL in keys(dst)
        dst[δL] .+= src[δL]
    end
    dst
end

function densitymatrix_compute_add!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64},
                                     μ::Float64, T::Real, exec::Parallel.Executor, progressbar=nothing; kwargs...)
    Parallel.configure_blas!(exec; verbose=false)

    # Lock-protected scalar accumulator for the kinetic energy. Lock contention
    # is negligible: the per-k work (eigen + density-matrix outer product)
    # dominates the addition by orders of magnitude.
    energy_lock = Threads.SpinLock()
    total_energy = Ref(0.0)

    # Same pattern for the progress bar: a SpinLock around `next!` is fine
    # because the bar update is microseconds vs millisecond+ per-k cost.
    progress_lock = Threads.SpinLock()

    Parallel.kspace_reduce!(ρs, ks, exec) do local_ρ, _scratch, j, k
        e = densitymatrix_serial_add_kpoint!(local_ρ, H, k, kweights[j]; μ=μ, T=T, kwargs...)
        Base.@lock energy_lock total_energy[] += e
        if progressbar !== nothing
            Base.@lock progress_lock ProgressMeter.next!(progressbar)
        end
    end
    total_energy[]
end

# Single-k kernel: factored out so both the kspace_reduce! body and the
# legacy serial path can share it.
function densitymatrix_serial_add_kpoint!(ρs::AbstractHops, H, k, w::Float64;
                                           μ::Float64=0.0, T::Real=0.01, kwargs...)
    Hk = H(k)
    ϵs, U = Eigen.geteigen(Hk; kwargs...)
    band = real.(ϵs)
    fd = fermidirac.(band .- μ; T=T)
    # ρ_k = U * Diagonal(fd) * U'  →  add  w * fourierphase(-k, δL) * ρ_k
    ρ_k = U * Diagonal(fd) * U'
    for δL in keys(ρs)
        ρs[δL] .+= w * fourierphase(-k, δL) .* ρ_k
    end
    real(w * sum(band .* fd))
end

