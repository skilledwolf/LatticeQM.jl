using ProgressMeter

import LatticeQM.Structure: regulargrid
import LatticeQM.Eigen
import LatticeQM.Spectrum
import LatticeQM.Parallel

const PROGRESSBAR_LDOS_DEFAULTLABEL = "LDOS"::String
const PROGRESSBAR_DENSITY_DEFAULTLABEL = "Density"::String
const PROGRESSBAR_DOS_DEFAULTLABEL = "DOS"::String

# Helper: build a Progress bar honoring the package-wide defaults.
_make_dos_progress(nks, label, hidebar) =
    Progress(nks; dt=Spectrum.PROGRESSBAR_MINTIME, desc=label, enabled=!hidebar)

##########################################################################################
# Density
##########################################################################################

function density_at_k!(n::AbstractVector{Float64}, spectrum_k, μ::Float64)
    ϵs, U = spectrum_k
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        if ϵ <= μ
            n .+= abs2.(ψ)
        end
    end
    n
end

function density(H, ks::AbstractMatrix{Float64}, μ::Float64=0.0; format=:dense, kwargs...)
    spectrum(k) = Eigen.geteigen(H(k); format=format)
    n = zeros(Float64, size(hamiltonian(ks[:, 1]), 1))
    density!(n, spectrum, ks, μ; kwargs...)
end

function density!(n::AbstractVector{Float64}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0;
                   multimode::Symbol=:auto, executor::Union{Nothing,Parallel.Executor}=nothing,
                   progress_label=PROGRESSBAR_DENSITY_DEFAULTLABEL,
                   hidebar=!Spectrum.PROGRESSBAR_SHOWDEFAULT)
    n .= 0
    L = size(ks, 2)
    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    Parallel.configure_blas!(exec; verbose=false)
    progressbar = _make_dos_progress(L, progress_label, hidebar)

    Parallel.kspace_reduce!(n, ks, exec; progress=progressbar) do local_n, _scratch, _j, k
        density_at_k!(local_n, spectrum(k), μ)
    end
    finish!(progressbar)
    n ./= L
end

##########################################################################################
# LDOS
##########################################################################################

function ldos!(n::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix, ωs::AbstractVector; Γ::Real=0.1)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        for ω in ωs
            n[:] .+= imag(abs2.(ψ) ./ (ω + 1.0im * Γ - ϵ))
        end
    end
    n[:] ./= size(ωs)
end

function ldos!(n::AbstractVector, H, ks::AbstractMatrix, ωs::AbstractVector;
                Γ::Real=0.1,
                multimode::Symbol=:auto,
                executor::Union{Nothing,Parallel.Executor}=nothing,
                progress_label=PROGRESSBAR_LDOS_DEFAULTLABEL,
                hidebar=!Spectrum.PROGRESSBAR_SHOWDEFAULT,
                format=:dense,
                kwargs...)
    L = size(ks, 2)
    n .= 0
    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    Parallel.configure_blas!(exec; verbose=false)
    progressbar = _make_dos_progress(L, progress_label, hidebar)

    Parallel.kspace_reduce!(n, ks, exec;
        scratch_factory = () -> (Hcache=Spectrum.bloch_buffer(H, ks; format=format),),
        progress = progressbar) do local_n, scratch, _j, k
        Hk = Spectrum.bloch!(scratch.Hcache, H, k)
        ϵs, U = Eigen.geteigen!(Hk; format=format, kwargs...)
        ldos!(local_n, ϵs, U, ωs; Γ=Γ)
    end
    finish!(progressbar)
    n .*= -1 / (L * π)
    n
end

ldos(H, ks, frequency::Real; kwargs...) = ldos(H, ks, [frequency]; kwargs...)
function ldos(H, ks, frequencies::AbstractVector; format=:sparse, kwargs...)
    n = zeros(Float64, dim(H, ks))
    ldos!(n, H, ks, frequencies; format=format, kwargs...)
    n
end

##########################################################################################
# DOS
##########################################################################################

"""
    getdos(h, emin::Float64, emax::Float64, num=500; kwargs...)

Computes the density of states of operator h(k) on the entire Brillouine zone,
discretized on a grid with \$ k_{lin} \\times k_{lin} \$ points. 
and for the frequencies ωs=(ωmin, ω2, ..., ωmax). The paremter \$\\Gamma\$ is the energy broadening.

Accepts the same kwargs as getdos(h, ωs; klin, Γ, kwargs...).

Note: the current implementation only works for a two-dimensional Brillouine zone.
Might change in the future, but for now use dos(h, ks, ω; Γ) syntax if needed.
"""
getdos(h, emin::Float64, emax::Float64, num::Int=500; kwargs...) = (Ωs=LinRange(emin, emax, num); (Ωs, getdos(h, Ωs; kwargs...)));

"""
    getdos(h, ωs; klin, Γ, kwargs...)

Computes the density of states of operator h(k) on the entire Brillouine zone,
discretized on a grid with \$ k_{lin} \\times k_{lin} \$ points. 
and for the frequencies ωs=(ω1, ω2, ...). The paremter \$\\Gamma\$ is the energy broadening.

Accepts the same kwargs as dos(h, ks, ω).

Note: the current implementation only works for a two-dimensional Brillouine zone.
Might change in the future, but for now use dos(h, ks, ω; Γ) syntax if needed.
"""
getdos(h, ωs; kwargs...) = (DOS=zero(ωs); getdos!(DOS, h, ωs; kwargs...))
function getdos!(DOS, h, ωs::AbstractVector; klin::Int, kwargs...)
    ks = regulargrid(nk=klin^2)
    getdos!(DOS, h, Vector(ωs), ks; kwargs...)

    DOS
end


getdos_dense(h, args...; kwargs...) = getdos(h, args...; format=:dense, kwargs...)
getdos_sparse(h, args...; kwargs...) = getdos(h, args...; format=:sparse, kwargs...)


"""
    getdos(h, ks, ω; Γ, parallel=true, format=:auto)

Computes the density of states of operator h(k) using the points ks=(k1,k2,...)
and for the frequencies ω=(ω1, ω2, ...). The paremter \$\\Gamma\$ is the energy broadening.

Mode can be :distributed or :serial, format can be :auto, :sparse or :dense.
"""
getdos(h, ωs::AbstractVector, ks, args...; kwargs...) = (DOS=zero(ωs); getdos!(DOS, h, ωs, ks, args...; kwargs...))

function getdos!(DOS, h, frequencies::AbstractVector, ks, args...;
                  multimode::Symbol=:auto,
                  executor::Union{Nothing,Parallel.Executor}=nothing,
                  progress_label=PROGRESSBAR_DOS_DEFAULTLABEL,
                  hidebar=!Spectrum.PROGRESSBAR_SHOWDEFAULT,
                  kwargs...)
    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    dos_compute!(DOS, h, frequencies, ks, args..., exec;
                 progress_label=progress_label, hidebar=hidebar, kwargs...)
    DOS
end

dos!(DOS, energies::AbstractVector, ωs; kwargs...) = foreach(x->dos!(DOS, x, ωs; kwargs...), energies)
function dos!(DOS, energy::Number, ωs; broadening::Number, weight::Number=1.0)
    DOS .+= imag.( weight.*1.0./(ωs .- 1.0im .* broadening .- energy) )
    DOS
end

import LatticeQM.Structure: Mesh, meshweights

# Single unified DOS kernel: chooses executor, runs the per-k accumulation,
# normalises by k-point count (for unweighted) or leaves weights as-is.
#
# `dos!` adds δ(ω - ϵ_k) (broadened) for each k to the local accumulator.
# Per-task accumulators are merged at the end by Parallel.kspace_reduce!.
#
# `kweights === nothing` is the unweighted path: per-k contributions are
# summed and divided by L at the end. Otherwise the weights are applied
# in-line and no post-normalisation is done (caller-supplied weights are
# expected to encode the BZ measure).
function dos_compute!(DOS, h, frequencies::AbstractVector,
                      ks::AbstractMatrix{<:Real}, exec::Parallel.Executor;
                      Γ::Number,
                      kweights::Union{Nothing,AbstractVector}=nothing,
                      progress_label=PROGRESSBAR_DOS_DEFAULTLABEL,
                      hidebar=!Spectrum.PROGRESSBAR_SHOWDEFAULT,
                      format=:dense, kwargs...)
    L = size(ks, 2)
    Parallel.configure_blas!(exec; verbose=false)
    progressbar = _make_dos_progress(L, progress_label, hidebar)
    Parallel.kspace_reduce!(DOS, ks, exec;
        scratch_factory = () -> (Hcache=Spectrum.bloch_buffer(h, ks; format=format),),
        progress = progressbar) do local_dos, scratch, j, k
        Hk = Spectrum.bloch!(scratch.Hcache, h, k)
        ϵs = Eigen.geteigvals!(Hk; format=format, kwargs...)
        w = kweights === nothing ? 1.0 : kweights[j]
        dos!(local_dos, ϵs, frequencies; broadening=Γ, weight=w)
    end
    finish!(progressbar)
    kweights === nothing && (DOS ./= L)
    DOS
end

dos_compute!(DOS, h, frequencies::AbstractVector,
             ks::AbstractMatrix, kweights::AbstractVector,
             exec::Parallel.Executor; kwargs...) =
    dos_compute!(DOS, h, frequencies, ks, exec; kweights=kweights, kwargs...)

dos_compute!(DOS, h, frequencies::AbstractVector, kgrid::Mesh,
             exec::Parallel.Executor; kwargs...) =
    dos_compute!(DOS, h, frequencies, kgrid.points, meshweights(kgrid), exec; kwargs...)

# using HCubature
#
# function dos_quad(ϵs::Function, energies::AbstractVector{Float64};  Γ::Float64, kwargs...)
#     fdim = length(energies)
#     f(k) = sum( imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) ) for ϵ=ϵs(k) )
#     xmin = [ 0.0, 0.0 ]
#     xmax = [ 1.0, 1.0 ]
#
#     (DOS, ERROR) = hcubature(f, xmin, xmax; kwargs...)
#
#     DOS/π, ERROR
# end

# using Cubature

# function dos_quad_parallel(ϵs::Function, energies::AbstractVector{T};  Γ::T, kwargs...) where {T<:Number}
#     fdim = length(energies)
#     f0(k) = sum( imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) ) for ϵ=ϵs(k) )
#     function f(ks, v)
#         v[:] = pmap(f0, eachcol(ks))
#     end
#     xmin = [ 0.0, 0.0 ]
#     xmax = [ 1.0, 1.0 ]

#     (DOS, ERROR) = hcubature_v(fdim, f, xmin, xmax; kwargs...)

#     DOS/π, ERROR
# end