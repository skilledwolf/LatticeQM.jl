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

# kspace_reduce! requires `output` to support `zero(output)` and accumulation.
# `Base.zero` for Hops produces a same-shape Hops with fresh zero matrices.
# The `_accumulate!` extension teaches `Parallel.kspace_reduce!` how to fold
# per-task partials back into the master output (Hops opts out of
# broadcasting via `Base.broadcastable(H::Hops) = Ref(H)`, so the default
# `output .+= partial` won't work for Hops).
Base.zero(h::Hops) = _densitymatrix_local(h)
function Parallel._accumulate!(out::AbstractHops, partial::AbstractHops)
    for δL in keys(out)
        out[δL] .+= partial[δL]
    end
    out
end

function densitymatrix_compute_add!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64},
                                     μ::Float64, T::Real, exec::Parallel.Executor, progressbar=nothing;
                                     format=:dense, kwargs...)
    Parallel.configure_blas!(exec; verbose=false)

    # Lock-protected scalar accumulator for the kinetic energy. Lock contention
    # is negligible: the per-k work (eigen + density-matrix outer product)
    # dominates the addition by orders of magnitude.
    energy_lock = Threads.SpinLock()
    total_energy = Ref(0.0)

    Parallel.kspace_reduce!(ρs, ks, exec;
        scratch_factory = () -> (Hcache=Spectrum.bloch_buffer(H, ks; format=format),),
        progress = progressbar) do local_ρ, scratch, j, k
        e = densitymatrix_kpoint!(local_ρ, scratch, H, k, kweights[j]; μ=μ, T=T, format=format, kwargs...)
        Base.@lock energy_lock total_energy[] += e
    end
    total_energy[]
end

# Single-k kernel: builds H(k) into the scratch buffer (zero-alloc for dense
# AbstractHops) and adds the contribution to the per-task accumulator.
function densitymatrix_kpoint!(ρs::AbstractHops, scratch, H, k, w::Float64;
                                μ::Float64=0.0, T::Real=0.01, format=:dense, kwargs...)
    Hk = Spectrum.bloch!(scratch.Hcache, H, k)
    ϵs, U = Eigen.geteigen!(Hk; format=format, kwargs...)
    band = real.(ϵs)
    fd = fermidirac.(band .- μ; T=T)
    ρ_k = U * Diagonal(fd) * U'
    for δL in keys(ρs)
        ρs[δL] .+= w * fourierphase(-k, δL) .* ρ_k
    end
    real(w * sum(band .* fd))
end

