# Convenient interfaces for hopping Hamiltonian

import ..Structure.Lattices: Lattice

"""
    opticalconductivityXX(frequencies, H, lat; kwargs...)

Convenience wrapper computing σ_xx(ω) via the Kubo formula. Builds current
operators from `H` and `lat` and evaluates at the provided `frequencies`.
"""
function opticalconductivityXX(frequencies::AbstractVector, H, lat::Lattice, args...; kwargs...)
    opticalconductivityXX(frequencies, 1, 1, H, lat, args...; kwargs...)
end

"""
    opticalconductivityXY(frequencies, H, lat; kwargs...)

Convenience wrapper computing σ_xy(ω). See `opticalconductivity(frequencies,
i, j, H, lat; ...)` for details.
"""
function opticalconductivityXY(frequencies::AbstractVector, H, lat::Lattice, args...; kwargs...)
    opticalconductivityXY(frequencies, 1, 2, H, lat, args...; kwargs...)
end

import ..Operators: getcurrentoperators

"""
    opticalconductivity(frequencies, i, j, H, lat; kwargs...)

Compute σ_ij(ω) for Hamiltonian `H(k)` using the Cartesian components `i, j` of
the current operator defined on `lat`. Internally obtains
`J = getcurrentoperators(lat, H)` and forwards to the low‑level routine.

Common keywords: `μ` (chemical potential), `Γ` (broadening), `T` (temperature),
plus diagonalization options forwarded to the eigen solver.
"""
function opticalconductivity(frequencies::AbstractVector, i::Int, j::Int, H, lat::Lattice, args...; kwargs...)
    J = getcurrentoperators(lat, H)
    opticalconductivity(frequencies, H, J[i], J[j], args...; kwargs...)
end

# More abstract interfaces

import ..Structure

"""
    opticalconductivity(frequencies, H, J1, J2; klin, kwargs...)

Compute σ(ω) on a regular `klin × klin` k‑grid using currents `J1, J2`.
Convenience front‑end to the explicit‑grid method.
"""
function opticalconductivity(frequencies::AbstractVector, H, J1, J2; klin, kwargs...)
    kgrid = Structure.regulargrid(;nk=klin^2)
    opticalconductivity(frequencies, H, J1, J2, kgrid; kwargs...)
end


import LatticeQM.Eigen
using SharedArrays
using ProgressMeter

import LinearAlgebra

"""
    opticalconductivity(frequencies, H, J1, J2, ks; μ=0.0, Γ=0.025, T=0.1, ...)

Evaluate the Kubo formula for the optical conductivity tensor at the set of
`frequencies` and k‑points `ks` (columns). `H`, `J1`, and `J2` are callable with
a k‑vector and return the Hamiltonian and current operators, respectively.

Returns a complex vector `σ(ω)` normalized by the number of k‑points. `μ` is the
chemical potential, `Γ` the phenomenological broadening and `T` the temperature
(all in the same energy units as `H`). Additional keyword arguments are
forwarded to the eigenvalue solver.
"""
function opticalconductivity(frequencies::AbstractVector, H, J1, J2, ks::AbstractMatrix; μ::Float64=0.0, Γ::Float64=0.025, T::Float64=0.1, kwargs...)

#     ks = points(ks)
    N = size(ks, 2)
    function spectrumf(k)
        Eigen.geteigen(H(k); kwargs...)
    end

#     OC = SharedArray(complex(zero(frequencies)))
    OC = complex(zero(frequencies))

#     @sync @showprogress 1 "Computing spectrum... " @distributed for j_=1:N
    @showprogress 1 "Computing spectrum... " for j_=1:N
        k = ks[:,j_]
        ϵs, U = spectrumf(k)

        kubo!(OC, frequencies, real(ϵs), U, J1(k), J2(k)'; μ=μ, Γ=Γ, T=T)
    end
    OC .= -1im .* OC ./ N

    OC
end
