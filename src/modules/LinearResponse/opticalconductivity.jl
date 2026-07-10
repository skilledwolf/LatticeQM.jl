# Convenient interfaces for hopping Hamiltonian

import ..Structure.Lattices: Lattice

"""
    opticalconductivityXX(frequencies, H, lat; kwargs...)

Convenience wrapper computing the xx component of the Kubo current-current
correlator χ_xx(ω). Builds current operators from `H` and `lat` and evaluates
at the provided `frequencies`.

Note: despite the name, this returns the raw correlator χ(ω), **not** the
optical conductivity σ(ω). See the low-level
`opticalconductivity(frequencies, H, J1, J2, ks; ...)` docstring for the
formula and the χ → σ conversion.
"""
function opticalconductivityXX(frequencies::AbstractVector, H, lat::Lattice, args...; kwargs...)
    opticalconductivityXX(frequencies, 1, 1, H, lat, args...; kwargs...)
end

"""
    opticalconductivityXY(frequencies, H, lat; kwargs...)

Convenience wrapper computing the xy component of the Kubo current-current
correlator χ_xy(ω) — not σ_xy(ω) itself; see the low-level
`opticalconductivity(frequencies, H, J1, J2, ks; ...)` docstring for the
χ → σ conversion. See `opticalconductivity(frequencies, i, j, H, lat; ...)`
for details.
"""
function opticalconductivityXY(frequencies::AbstractVector, H, lat::Lattice, args...; kwargs...)
    opticalconductivityXY(frequencies, 1, 2, H, lat, args...; kwargs...)
end

import ..Operators: getcurrentoperators

"""
    opticalconductivity(frequencies, i, j, H, lat; kwargs...)

Compute the Kubo current-current correlator χ_ij(ω) for Hamiltonian `H(k)`
using the Cartesian components `i, j` of the current operator defined on
`lat`. Internally obtains `J = getcurrentoperators(lat, H)` and forwards to
the low‑level routine, whose docstring gives the formula and the conversion
from the returned χ(ω) to the optical conductivity σ(ω).

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

Compute the Kubo current-current correlator χ(ω) on a regular `klin × klin`
k‑grid using currents `J1, J2`. Convenience front‑end to the explicit‑grid
method, whose docstring gives the formula and the χ(ω) → σ(ω) conversion.
"""
function opticalconductivity(frequencies::AbstractVector, H, J1, J2; klin, kwargs...)
    kgrid = Structure.regulargrid(;nk=klin^2)
    opticalconductivity(frequencies, H, J1, J2, kgrid; kwargs...)
end


import LatticeQM.Eigen
import LatticeQM.Parallel
import LatticeQM.Spectrum
using ProgressMeter

import LinearAlgebra

"""
    opticalconductivity(frequencies, H, J1, J2, ks; μ=0.0, Γ=0.025, T=0.1, ...)

Evaluate the Kubo current-current correlator at the set of `frequencies` and
k‑points `ks` (columns). `H`, `J1`, and `J2` are callable with a k‑vector and
return the Hamiltonian and current operators, respectively.

Returns a complex vector `χ(ω)` normalized by the number of k‑points `N`:

    χ(ω) = -(i/N) Σ_k Σ_{m≠n} (f_m - f_n) ⟨m|J1|n⟩⟨n|J2†|m⟩ / (ϵ_m - ϵ_n - ω - iΓ)

where `f = fermidirac(ϵ; μ, T)` and `|m⟩, ϵ_m` are the Bloch eigenstates at k.
This is **not** yet the optical conductivity. Obtain σ(ω) by subtracting the
static (diamagnetic) piece and dividing by the photon energy:

    σ(ω) = -(χ(ω) - χ(0)) / (ω + iΓ)

as done in `test/test_linearresponse.jl` and
`extra/examples/graphene/opticalconductivity.jl`.

`μ` is the chemical potential, `Γ` the phenomenological broadening and `T` the
temperature (all in the same energy units as `H`). Additional keyword
arguments are forwarded to the eigenvalue solver.
"""
function opticalconductivity(frequencies::AbstractVector, H, J1, J2, ks::AbstractMatrix;
                              μ::Float64=0.0, Γ::Float64=0.025, T::Float64=0.1,
                              multimode::Symbol=:auto,
                              executor::Union{Nothing,Parallel.Executor}=nothing,
                              progress_label="Optical conductivity",
                              hidebar=false,
                              format=:dense,
                              kwargs...)
    N = size(ks, 2)
    OC = complex(zero(frequencies))

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    Parallel.configure_blas!(exec; verbose=false)

    progressbar = ProgressMeter.Progress(N; dt=1, desc=progress_label, enabled=!hidebar)

    Parallel.kspace_reduce!(OC, ks, exec;
        scratch_factory = () -> (Hcache=Spectrum.bloch_buffer(H, ks; format=format),),
        progress = progressbar) do local_OC, scratch, _j, k
        Hk = Spectrum.bloch!(scratch.Hcache, H, k)
        ϵs, U = Eigen.geteigen!(Hk; format=format, kwargs...)
        kubo!(local_OC, frequencies, real(ϵs), U, J1(k), J2(k)'; μ=μ, Γ=Γ, T=T)
    end
    ProgressMeter.finish!(progressbar)
    OC .*= -1im / N
    OC
end
