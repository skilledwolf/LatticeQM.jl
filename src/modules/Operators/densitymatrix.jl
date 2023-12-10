using Distributed
using SharedArrays
using ProgressMeter
# import Tullio # for tensor contractions
# import TensorOperations # for tensor contractions

import ..Utils: fermidirac
import ..TightBinding: Hops, AbstractHops
import ..Spectrum: dim
import ..TightBinding: fourierphase
import ..Structure: Mesh, meshweights

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
        densitymatrix_distributed!(ρs, H, ks, kweights, μ; kwargs...) 
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

function densitymatrix_distributed!(ρs::AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, progressmin::Int=20, kwargs...)
    L = size(ks, 2)

    energies = SharedArray(zeros(Float64, L))
    δLs = collect(keys(ρs))

    M = length(δLs); Ns = size(H)
    ρ0 = SharedArray(zeros(ComplexF64, Ns..., M))

    fourierphases = [fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]  

    @sync @showprogress dt = Spectrum.PROGRESSBAR_MINTIME desc = PROGRESSBAR_DENSMAT_DEFAULTLABEL enabled = Spectrum.PROGRESSBAR_SHOWDEFAULT @distributed for i_ = 1:L
        k = ks[:, i_]
        ϵs_k, U_k = Eigen.geteigen(H(k); kwargs...)
        fd_k = fermidirac.(real.(ϵs_k .- μ); T=T)

        w = kweights[i_]
        phases = fourierphases[:, i_] # [fourierphase(-k, δL) for δL in δLs]

        U_k .= U_k * Diagonal(fd_k) * U_k' # not sure if this works correctly in-place?

        for j_ = axes(ρ0, 3)
            ρ0[:, :, j_] .+= w .* phases[j_] .* U_k
        end

        energies[i_] = real(w * sum(ϵs_k .* fd_k))
        # GC.gc() # force garbage collection
    end

    flexibleformat!(ρs, ρ0, δLs)
    sum(energies) # return the kinetic part of the gs energy
end


import ..Spectrum

function densitymatrix_serial!(ρs::TightBinding.AbstractHops, H, ks::AbstractMatrix{Float64}, kweights::AbstractVector{Float64}, μ::Float64=0.0; T::Real=0.01, kwargs...)

    L = size(ks, 2)
    energies = zeros(Float64, L)

    δLs = keys(ρs)
    for δL in δLs
        ρs[δL] .= 0.0
    end

    fourierphases = [TightBinding.fourierphase(-k, δL) for δL in δLs, k in eachcol(ks)]
    ρ_k = zeros(ComplexF64, size(H))

    @showprogress dt=Spectrum.PROGRESSBAR_MINTIME desc=PROGRESSBAR_DENSMAT_DEFAULTLABEL enabled=Spectrum.PROGRESSBAR_SHOWDEFAULT for i_ = 1:L
        k = ks[:, i_]
        w_k = kweights[i_]

        ϵs_k, U_k = Eigen.geteigen(H(k); kwargs...)
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