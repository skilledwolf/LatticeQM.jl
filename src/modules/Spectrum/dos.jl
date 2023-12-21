using Distributed
# using SharedArrays
import LatticeQM.Structure: regulargrid

using Distributed
using ProgressMeter

import LatticeQM.Eigen
import LatticeQM.Spectrum

const PROGRESSBAR_LDOS_DEFAULTLABEL = "LDOS"::String

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

function density!(n::AbstractVector{Float64}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)
    n .= zero(n)
    L = size(ks, 2)

    n .= @distributed (+) for j = 1:L # @todo: this should be paralellized
        n0 = zero(n)    ## <-- it annoys me that I don't know how to get around this allocation
        density_at_k!(n0, spectrum(ks[:, j]), μ)
        n0 .= n0 ./ L
    end
    n
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

function ldos!(n::AbstractVector, H, ks::AbstractMatrix, ωs::AbstractVector; Γ::Real=0.1, progress_label=PROGRESSBAR_LDOS_DEFAULTLABEL, kwargs...)
    L = size(ks, 2)


    n[:] = @sync @showprogress dt = Spectrum.PROGRESSBAR_MINTIME desc = progress_label enabled = Spectrum.PROGRESSBAR_SHOWDEFAULT @distributed (+) for j = 1:L
        n0 = zero(n)
        ϵs, U = Eigen.geteigen(H(ks[:, j]); kwargs...)
        ldos!(n0, ϵs, U, ωs; Γ=Γ)
        n0
    end
    n[:] .= -n ./ L ./ π
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
function getdos!(DOS, h, frequencies::AbstractVector, ks, args...; parallel=true, mode=:distributed, kwargs...)
    if nprocs()<2 || mode!=:distributed
        parallel=false
    end
    
    if parallel
        DOS = dos_parallel!(DOS, h, frequencies, ks, args...; kwargs...)
    else
        DOS = dos_serial!(DOS, h, frequencies, ks, args...; kwargs...)
    end

    DOS
end

dos!(DOS, energies::AbstractVector, ωs; kwargs...) = foreach(x->dos!(DOS, x, ωs; kwargs...), energies)
function dos!(DOS, energy::Number, ωs; broadening::Number, weight::Number=1.0)
    DOS .+= imag.( weight.*1.0./(ωs .- 1.0im .* broadening .- energy) )
    DOS
end

function dos_serial!(DOS, h, frequencies::AbstractVector, ks::AbstractMatrix{<:Real}; Γ::Number, kwargs...)
    L = size(ks,2)
    function ϵs(k)
        let h0 = h(k), kwargs = kwargs
            Eigen.geteigvals(h0; kwargs...)
        end
    end

    @showprogress 6 "Computing DOS... " for k=eachcol(ks) # j=1:L
        dos!(DOS, ϵs(k), frequencies; broadening=Γ)
    end
    DOS ./= L
    DOS
end

function dos_parallel!(DOS, h, frequencies::AbstractVector, ks::AbstractMatrix; Γ::Number, kwargs...)
    L = size(ks,2)
    function ϵs(k)
        let h0 = h(k), kwargs = kwargs
            Eigen.geteigvals(h0; kwargs...)
        end
    end


    DOS0 = @sync @showprogress 6 "Computing DOS... " @distributed (+) for j=1:L # over ks
        tmp = zero(DOS)
        dos!(tmp, ϵs(ks[:,j]), frequencies; broadening=Γ)

        tmp
    end
    DOS[:] += (DOS0 / L)[:]

    DOS
end

function dos_serial!(DOS, h, frequencies::AbstractVector, ks::AbstractMatrix, kweights::AbstractVector; Γ::Number, kwargs...)
    L = size(ks,2)
    function ϵs(k)
        let h0=h(k), kwargs=kwargs
            Eigen.geteigvals(h0; kwargs...)
        end
    end

    @showprogress 6 "Computing DOS... " for (k,w)=zip(eachcol(ks),kweights) # j=1:L
        dos!(DOS, ϵs(k), frequencies; broadening=Γ, weight=w)
    end
    # DOS ./= L
    DOS
end

function dos_parallel!(DOS, h, frequencies::AbstractVector, ks::AbstractMatrix, kweights::AbstractVector; Γ::Number, kwargs...)
    L = size(ks,2)
    function ϵs(k)
        let h0 = h(k), kwargs = kwargs
            Eigen.geteigvals(h0; kwargs...)
        end
    end


    DOS0 = @sync @showprogress 6 "Computing DOS... " @distributed (+) for j=1:L # over ks
        tmp = zero(DOS)
        dos!(tmp, ϵs(ks[:,j]), frequencies; broadening=Γ, weight=kweights[j])

        tmp
    end
    # DOS[:] += (DOS0 / L)[:]
    DOS[:] += DOS0[:]

    DOS
end

import LatticeQM.Structure: Mesh, meshweights

function dos_serial!(DOS, h, frequencies::AbstractVector, kgrid::Mesh; kwargs...)
    ks = kgrid.points
    kweights = meshweights(kgrid)
    dos_serial!(DOS, h, frequencies, ks, kweights; kwargs...)
end

function dos_parallel!(DOS, h, frequencies::AbstractVector, kgrid::Mesh; kwargs...)
    ks = kgrid.points
    kweights = meshweights(kgrid)
    dos_parallel!(DOS, h, frequencies, ks, kweights; kwargs...)
end

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