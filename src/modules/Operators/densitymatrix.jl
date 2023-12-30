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

densitymatrix(ϵ::Number, ψ::AbstractVector; T::Real = 0.01) = fermidirac(real(ϵ); T = T) .* transpose(ψ * ψ')

function densitymatrix!(ρ_k::AbstractMatrix, fd_k::AbstractVector, U_k::AbstractMatrix)
    mul!(ρ_k, U_k, Diagonal(fd_k) * U_k')
    ρ_k
end

function getdensitymatrices_Ls(phases::AbstractVector, fd::AbstractVector, U::AbstractMatrix, w::Number=1.0)
    D, N = size(U)
    mat = Array{ComplexF64}(undef, (D, D))
    densitymatrix!(mat, fd, U)
    ρ0 = Array{ComplexF64}(undef, (D, D, size(phases,1)))
    for j_ = axes(phases,1)
        @views ρ0[:, :, j_] .= w .* phases[j_] .* mat
    end
    ρ0
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

function getdensitymatrix!(ρs::Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; multimode=:serial, kwargs...)

    if multimode == :distributed && nprocs() > 1
        # densitymatrix_distributed!(ρs, H, ks, kweights, μ; kwargs...) 
        # densitymatrix_distributed_views!(ρs, H, ks, kweights, μ; kwargs...)
        densitymatrix_pmap_views!(ρs, H, ks, kweights, μ; kwargs...)
        # densitymatrix_dagger!(ρs, H, ks, kweights, μ; kwargs...) 
    else
        densitymatrix_serial!(ρs, H, ks, kweights, μ; kwargs...)
    end
end

###################################################################################################
###################################################################################################
###################################################################################################

import LatticeQM.Spectrum
import LatticeQM.Eigen
import LatticeQM.TightBinding: efficientzero, flexibleformat!, fourierphase

# @todo: phase out this function, redefine the shared-hops type to use views in the first place
#        and get rid of the ugly hack below
# function densitymatrix_distributed!(ρs::Hops, H::T1, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...) where {T1}
#     L = size(ks, 2)

#     energies = SharedArray(zeros(Float64, L))
#     δLs = collect(keys(ρs))

#     M = length(δLs); Ns = size(H)
#     ρ0 = SharedArray(zeros(ComplexF64, Ns..., M))

#     fourierphases = [fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]  

#     @sync @showprogress dt = Spectrum.PROGRESSBAR_MINTIME desc = PROGRESSBAR_DENSMAT_DEFAULTLABEL enabled = Spectrum.PROGRESSBAR_SHOWDEFAULT @distributed for i_ = 1:L
        
#         # h0::Matrix{ComplexF64} = H(ks[:, i_])
#         # ϵs_k, U_k = Eigen.geteigen(h0; kwargs...)
#         ϵs_k, U_k = Eigen.geteigen!(H(ks[:, i_]); kwargs...)
        
#         fd_k = fermidirac.(real.(ϵs_k .- μ); T=T)

#         w = kweights[i_]
#         phases = fourierphases[:, i_] 

#         U_k .= U_k * Diagonal(fd_k) * U_k' 

#         for j_ = axes(ρ0, 3)
#             ρ0[:, :, j_] .+= w .* phases[j_] .* U_k
#         end
#         U_k = nothing

#         energies[i_] = real(w * sum(ϵs_k .* fd_k))
#     end

#     flexibleformat!(ρs, ρ0, δLs)
#     ρ0 = nothing
#     sum(energies) # return the kinetic part of the gs energy
# end

import LatticeQM.TightBinding
import LatticeQM.TightBinding: SharedDenseHops, Hops, SubarrayHops

function densitymatrix_pmap_views!(ρs::SharedDenseHops, args...; kwargs...)
    # ρs_views = TightBinding.Hops(Dict(L => view(M, :, :) for (L, M) in ρs))
    ρs_views = TightBinding.gethopsview(ρs)
    densitymatrix_pmap_views!(ρs_views, args...; kwargs...)
end

function densitymatrix_pmap_views!(ρs::SubarrayHops, H::T1, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...) where {T1}
    L = size(ks, 2)

    δLs = collect(keys(ρs))

    for δL in δLs
        ρs[δL] .= 0.0
    end

    fourierphases = [fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]

    kinetic_energies = @showprogress dt = Spectrum.PROGRESSBAR_MINTIME desc = PROGRESSBAR_DENSMAT_DEFAULTLABEL enabled = Spectrum.PROGRESSBAR_SHOWDEFAULT pmap(1:L) do i_
        ϵs_k, U_k = Eigen.geteigen!(H(ks[:, i_]); kwargs...)
        fd_k = fermidirac.(real.(ϵs_k .- μ); T=T)

        w = kweights[i_]
        phases = fourierphases[:, i_] # [fourierphase(-k, δL) for δL in δLs]

        U_k .= U_k * Diagonal(fd_k) * U_k' # not sure if this works correctly in-place?

        for (n_, δL) in enumerate(δLs)
            # CRUCIAL NOTE: for sharedarrays in distributed mode, this will 
            # only work as expected if we reference the hops through views!!!
            # Hops cannot directly contain sharedarrays (at least in julia 1.10)
            @. ρs[δL] += w * phases[n_] * U_k
        end

        real(w * sum(ϵs_k .* fd_k))
    end

    sum(kinetic_energies) # return the kinetic part of the gs energy
end

function densitymatrix_distributed_views!(ρs::SharedDenseHops, args...; kwargs...)
    # ρs_views = TightBinding.Hops(Dict(L => view(M, :, :) for (L, M) in ρs))
    ρs_views = TightBinding.gethopsview(ρs)
    densitymatrix_distributed_views!(ρs_views, args...; kwargs...)
end 

function densitymatrix_distributed_views!(ρs::SubarrayHops, H::T1, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...) where {T1}
    L = size(ks, 2)
    δLs = collect(keys(ρs))

    for δL in δLs
        ρs[δL] .= 0.0
    end

    fourierphases = [fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]

    kinetic_energy = @showprogress dt = Spectrum.PROGRESSBAR_MINTIME desc = PROGRESSBAR_DENSMAT_DEFAULTLABEL enabled = Spectrum.PROGRESSBAR_SHOWDEFAULT @distributed (+) for i_ = 1:L

        ϵs_k, U_k = Eigen.geteigen!(H(ks[:, i_]); kwargs...)
        fd_k = fermidirac.(real.(ϵs_k .- μ); T=T)

        w = kweights[i_]
        phases = fourierphases[:, i_] # [fourierphase(-k, δL) for δL in δLs]

        U_k .= U_k * Diagonal(fd_k) * U_k' # not sure if this works correctly in-place?

        for (n_, δL) in enumerate(δLs)
            # CRUCIAL NOTE: for sharedarrays in distributed mode, this will 
            # only work as expected if we reference the hops through views!!!
            # Hops cannot directly contain sharedarrays (at least in julia 1.10)
            @. ρs[δL] += w * phases[n_] * U_k
        end
        U_k = nothing

        real(w * sum(ϵs_k .* fd_k))
    end

    kinetic_energy # return the kinetic part of the gs energy
end


import LatticeQM.Spectrum

function densitymatrix_serial!(ρs::TightBinding.Hops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, kwargs...)

    L = size(ks, 2)
    energies = zeros(Float64, L)

    δLs = keys(ρs)
    for δL in δLs
        ρs[δL] .= 0.0
    end

    fourierphases = [TightBinding.fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]
    ρ_k = zeros(ComplexF64, size(H))

    @showprogress dt=Spectrum.PROGRESSBAR_MINTIME desc=PROGRESSBAR_DENSMAT_DEFAULTLABEL enabled=Spectrum.PROGRESSBAR_SHOWDEFAULT for i_ = 1:L
        w_k = kweights[i_]

        ϵs_k, U_k = Eigen.geteigen!(H(ks[:, i_]); kwargs...)
        fd_k = fermidirac.(real.(ϵs_k .- μ); T=T)
        phases = fourierphases[:, i_]

        mul!(ρ_k, U_k, Diagonal(fd_k) * U_k')

        for (n_, δL) in enumerate(δLs)
            ρs[δL] .+= w_k .* ρ_k .* phases[n_]
        end

        energies[i_] = real(w_k * sum(ϵs_k .* fd_k))
    end

    sum(energies) # return the kinetic part of the gs energy
end


# import Dagger
# function densitymatrix_dagger!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
#     L = size(ks, 2)
#     energies = Dagger.@mutable zeros(Float64, L)
#     δLs = collect(keys(ρs))

#     M = length(δLs)
#     Ns = size(H)
#     ρsMat = Dagger.@mutable zeros(ComplexF64, Ns..., M)

#     function getspectrum(i_)
#         k = ks[:, i_]
#         Eigen.geteigen(H(k); kwargs...)
#     end

#     # Launch spectrum calculations
#     spectrum_tasks = []
#     for i_ = 1:L
#         push!(spectrum_tasks, (i_, Dagger.@spawn getspectrum(i_)))
#     end

#     # Launch addition tasks
#     function addToRho!(ρsMat, energies, spectrum_k)
#         (i_, (ϵs, U)) = fetch(spectrum_k)
#         w = kweights[i_]
#         k = ks[:, i_]

#         fd = fermidirac.(real.(ϵs .- μ); T=T)
#         phases = [fourierphase(-k, δL) for δL in δLs]

#         ρ_k = U * Diagonal(fd) * U'
#         for j_ = 1:M
#             ρsMat[:, :, j_] .+= w .* phases[j_] .* ρ_k
#         end

#         energies[i_] = real(w * sum(ϵs .* fd))
#         nothing
#     end

#     addition_taks = [Dagger.@spawn addToRho!(ρsMat, energies, spectrum_k) for spectrum_k in spectrum_tasks]
#     wait.(addition_taks)

#     flexibleformat!(ρs, collect(ρsMat), δLs)
#     sum(collect(energies)) # return the kinetic part of the gs energy
# end