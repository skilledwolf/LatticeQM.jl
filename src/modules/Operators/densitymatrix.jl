using Distributed
using SharedArrays
using ProgressMeter
# import Tullio # for tensor contractions
# import TensorOperations # for tensor contractions

import LatticeQM.Utils: fermidirac
import LatticeQM.TightBinding: Hops, AbstractHops
import LatticeQM.Spectrum: dim
import LatticeQM.TightBinding: fourierphase
import LatticeQM.Structure: Mesh, meshweights

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

function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; multimode=:serial, progress_label=PROGRESSBAR_DENSMAT_DEFAULTLABEL, hidebar=!Spectrum.PROGRESSBAR_SHOWDEFAULT, kwargs...)

    for δL in keys(ρs)
        ρs[δL] .= 0
    end

    progressbar = Progress(size(ks, 2); dt=Spectrum.PROGRESSBAR_MINTIME, desc=progress_label, enabled=!hidebar)

    out::Float64 = 0
    if multimode == :distributed && nprocs() > 1
        # out = densitymatrix_serial_add!(ρs, H, ks, kweights, progressbar; μ=μ, kwargs...)
        out = densitymatrix_distributed_add!(ρs, H, ks, kweights, progressbar; μ=μ, kwargs...)
        # out = densitymatrix_pmap_add!(ρs, H, ks, kweights, μ, progressbar; kwargs...)
    else
        if multimode == :multithreaded
            @warn "Multithreaded mode is not yet supported/implemented. Falling back to serial mode."
        end
        out = densitymatrix_serial_add!(ρs, H, ks, kweights, progressbar; μ=μ, kwargs...)
    end
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
import LatticeQM.TightBinding: SharedDenseHops, Hops, SubarrayHops


densitymatrix_distributed_add!(ρs::SharedDenseHops, args...; kwargs...) = densitymatrix_distributed_add!(TightBinding.gethopsview(ρs), args...; kwargs...)

function _densitymatrix_allocate_local(δLs, sizes)
    isempty(δLs) && return TightBinding.Hops(Dict{Vector{Int}, Matrix{ComplexF64}}())
    key_type = eltype(δLs)
    data = Dict{key_type, Matrix{ComplexF64}}()
    for (δL, sz) in zip(δLs, sizes)
        data[δL] = zeros(ComplexF64, sz)
    end
    TightBinding.Hops(data)
end

function _densitymatrix_compute_chunk(H, ks, kweights, chunk, δLs, sizes, μ, T, kwargs::NamedTuple)
    local_hops = _densitymatrix_allocate_local(δLs, sizes)
    if isempty(δLs)
        return local_hops, 0.0, length(chunk)
    end
    @views energy = densitymatrix_serial_add!(local_hops, H, ks[:, chunk], kweights[chunk]; μ=μ, T=T, kwargs...)
    local_hops, energy, length(chunk)
end

function _densitymatrix_distributed_add!(ρs::SubarrayHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, progressbar;
        μ::Float64=0.0, T::Real = 0.01, num_chunks=:auto, kwargs...)

    kwargs_nt = (; kwargs...)
    workers_list = workers()

    if isempty(workers_list)
        return densitymatrix_serial_add!(ρs, H, ks, kweights, progressbar; μ=μ, T=T, kwargs_nt...)
    end

    nks = size(ks, 2)
    num_chunks = num_chunks == :auto ? length(workers_list) : num_chunks
    num_chunks = max(num_chunks, 1)
    chunks = Spectrum.compute_chunks(nks, num_chunks)

    δLs = collect(keys(ρs))
    sizes = [size(ρs[δL]) for δL in δLs]

    futures = Vector{Future}(undef, length(chunks))
    for (i, chunk) in enumerate(chunks)
        worker_id = workers_list[mod1(i, length(workers_list))]
        futures[i] = @spawnat worker_id _densitymatrix_compute_chunk(H, ks, kweights, chunk, δLs, sizes, μ, T, kwargs_nt)
    end

    total_energy = 0.0
    for fut in futures
        local_hops, energy, chunk_len = fetch(fut)
        total_energy += energy
        for δL in δLs
            @. ρs[δL] += local_hops[δL]
        end
        if progressbar !== nothing
            for _ in 1:chunk_len
                ProgressMeter.next!(progressbar)
            end
        end
    end

    total_energy
end

function densitymatrix_distributed_add!(ρs::SubarrayHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}; kwargs...)
    _densitymatrix_distributed_add!(ρs, H, ks, kweights, nothing; kwargs...)
end

function densitymatrix_distributed_add!(ρs::SubarrayHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, progressbar; kwargs...)
    _densitymatrix_distributed_add!(ρs, H, ks, kweights, progressbar; kwargs...)
end

function densitymatrix_serial_add!(ρs::TightBinding.Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, progressbar=nothing; T::Real=0.01, μ::Float64=0.0, progress_channel = nothing, kwargs...)

    L = size(ks, 2)
    δLs = keys(ρs)
    energies = Vector{Float64}(undef, L)

    fd_cache = Vector{Float64}(undef, size(H, 1))
    band_cache = Vector{ComplexF64}(undef, size(H, 1))

    fourierphases = [TightBinding.fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]
    Hcache = Matrix{ComplexF64}(undef, size(H))
    Ucache = Matrix{ComplexF64}(undef, size(H))

    for i_ = axes(ks, 2)
        w = kweights[i_]
        Hcache .= H(ks[:, i_])
        band_cache, Ucache = Eigen.geteigen!(Hcache; kwargs...)
        fd_cache .= fermidirac.(real.(band_cache .- μ); T=T)

        # mul!(ρ_k, U_k, Diagonal(fd_cache) * U_k')
        Ucache .= Ucache * Diagonal(fd_cache) * Ucache'

        for (n_, δL) in enumerate(δLs)
            @. ρs[δL] += w * fourierphases[n_, i_] * Ucache #ρ_k
        end

        energies[i_] = real(w * sum(band_cache .* fd_cache))

        !isnothing(progress_channel) && put!(progress_channel, true)
        !isnothing(progressbar) && ProgressMeter.next!(progressbar)
    end
    Hcache = nothing
    Ucache = nothing
    band_cache = nothing
    fd_cache = nothing
    GC.gc() # eagerly collect garbage, especially for the local caches

    sum(energies) # return the kinetic part of the gs energy
end
