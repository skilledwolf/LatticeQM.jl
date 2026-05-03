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

# Per-task scratch for the density-matrix loop: dense H buffer + dense ρ_k
# buffer + temporary U·Diagonal(fd) buffer. Reused across every k handled by
# the same task — saves one ~N² ComplexF64 allocation per k, which on TBG-N=11
# (1588 orbitals, 16 k-points) is ~80 MB / SCF iteration of GC pressure.
#
# For sparse format with `num_bands < N`, U has shape (N, num_bands) and the
# UD intermediate matches; ρ_k is still N×N because `U Diagonal(fd) U'` is
# a rank-`num_bands` projector embedded in the full N×N space.
function _densitymatrix_scratch(ρs, H, ks; format=:dense, kwargs...)
    Hbuf = Spectrum.bloch_buffer(H, ks; format=format)
    N = Spectrum.dim(H, ks)
    nb = format === :sparse ? get(kwargs, :num_bands, N)::Int : N
    ρ_k = Matrix{ComplexF64}(undef, N, N)
    UD = Matrix{ComplexF64}(undef, N, nb)
    return (Hcache=Hbuf, ρ_k=ρ_k, UD=UD)
end

function densitymatrix_compute_add!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64},
                                     μ::Float64, T::Real, exec::Parallel.Executor, progressbar=nothing;
                                     format=:dense, kwargs...)
    Parallel.configure_blas!(exec; verbose=false)

    # `kspace_reduce!`'s `scalar_init` collects per-k kinetic-energy
    # contributions (the body's return value) and sums them across tasks /
    # workers. This sidesteps the closure-capture pitfall under distributed
    # mode (a `Ref` would be serialised per worker and the master would
    # never see the writes — see Parallel.kspace_reduce! docstring).
    Parallel.kspace_reduce!(ρs, ks, exec;
        scratch_factory = () -> _densitymatrix_scratch(ρs, H, ks; format=format, kwargs...),
        scalar_init = 0.0,
        progress = progressbar) do local_ρ, scratch, j, k
        densitymatrix_kpoint!(local_ρ, scratch, H, k, kweights[j]; μ=μ, T=T, format=format, kwargs...)
    end
end

# Single-k kernel: builds H(k) into the scratch buffer (zero-alloc for dense
# AbstractHops) and adds the contribution to the per-task accumulator.
# `scratch.UD` and `scratch.ρ_k` are preallocated dense buffers reused across
# k — replaces the temporary `U * Diagonal(fd) * U'` allocation that
# previously freshly allocated O(N²) per k.
function densitymatrix_kpoint!(ρs::AbstractHops, scratch, H, k, w::Float64;
                                μ::Float64=0.0, T::Real=0.01, format=:dense, kwargs...)
    Hk = Spectrum.bloch!(scratch.Hcache, H, k)
    ϵs, U = Eigen.geteigen!(Hk; format=format, kwargs...)
    band = real.(ϵs)
    fd = fermidirac.(band .- μ; T=T)
    # ρ_k = U * Diagonal(fd) * U' — done as two mul!s through a preallocated
    # intermediate UD = U * Diagonal(fd). For dense U with shape N×N this
    # eliminates two N²-sized temporaries per k.
    mul!(scratch.UD, U, Diagonal(fd))
    mul!(scratch.ρ_k, scratch.UD, U')
    for δL in keys(ρs)
        ρs[δL] .+= w * fourierphase(-k, δL) .* scratch.ρ_k
    end
    real(w * sum(band .* fd))
end

