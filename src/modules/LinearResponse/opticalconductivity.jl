# Convenient interfaces for hopping Hamiltonian

import ..Structure.Lattices: Lattice

function opticalconductivityXX(frequencies::AbstractVector, H, lat::Lattice, args...; kwargs...)
    opticalconductivityXX(frequencies, 1, 1, H, lat, args...; kwargs...)
end

function opticalconductivityXY(frequencies::AbstractVector, H, lat::Lattice, args...; kwargs...)
    opticalconductivityXY(frequencies, 1, 2, H, lat, args...; kwargs...)
end

import ..Operators: getcurrentoperators

function opticalconductivity(frequencies::AbstractVector, i::Int, j::Int, H, lat::Lattice, args...; kwargs...)
    J = getcurrentoperators(lat, H)
    opticalconductivity(frequencies, H, J[i], J[j], args...; kwargs...)
end

# More abstract interfaces

import ..Structure

function opticalconductivity(frequencies::AbstractVector, H, J1, J2; klin, kwargs...)
    kgrid = Structure.regulargrid(;nk=klin^2)
    opticalconductivity(frequencies, H, J1, J2, kgrid; kwargs...)
end


import ..Spectrum
using SharedArrays
using ProgressMeter

import LinearAlgebra

function opticalconductivity(frequencies::AbstractVector, H, J1, J2, ks::AbstractMatrix; μ::Float64=0.0, Γ::Float64=0.025, T::Float64=0.1, kwargs...)

#     ks = points(ks)
    N = size(ks, 2)
    function spectrumf(k)
        Spectrum.spectrum(H(k); kwargs...)
    end

#     OC = SharedArray(complex(zero(frequencies)))
    OC = complex(zero(frequencies))

#     @sync @showprogress 1 "Computing spectrum... " @distributed for j_=1:N
    @showprogress 1 "Computing spectrum... " for j_=1:N
        k = ks[:,j_]
        ϵs, U = spectrumf(k)

        kubo!(OC, frequencies, real(ϵs), U, J1(k), J2(k); μ=μ, Γ=Γ, T=T)
    end
    OC .= OC ./ N

#     (OC.-OC[1]) ./ (frequencies .+ 1im*Γ) # Calculate the optical conductivity tensor
    OC
end

