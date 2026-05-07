import LatticeQM.Operators
using LatticeQM.Operators: addchemicalpotential!

using LatticeQM.TightBinding: Hops, zerokey, hopdim
using SparseArrays#: spzeros

import LatticeQM.Spectrum
import LatticeQM.Utils
import LatticeQM.Structure


Spectrum.dim(h::BdGOperator, args...) = Spectrum.dim(h.h, args...)

function Operators.addchemicalpotential!(H::BdGOperator, μ::Real)
    d = hopdim(H)
    addhops!(H, BdGOperator(Hops(zerokey(H.h)=>complex(spzeros(d,d)+μ*I))))
end

import ..Structure: Lattices

function _addchemicalpotential_lattice_bdg!(H::BdGOperator, lat::Lattices.Lattice, μ; kwargs...)
    d = hopdim(H)
    chemhops = Hops(zerokey(H)=>zeros(ComplexF64,d,d))
    Operators.addchemicalpotential!(chemhops, lat, μ; kwargs...)

    # addhops!(H, BdGOperator(chemhops))
    addelectronsector!(H, chemhops)
end
Operators.addchemicalpotential!(H::BdGOperator, lat::Lattices.Lattice, μ::Float64; kwargs...) =
    _addchemicalpotential_lattice_bdg!(H, lat, μ; kwargs...)
Operators.addchemicalpotential!(H::BdGOperator, lat::Lattices.Lattice, μ::AbstractVector{<:Float64}; kwargs...) =
    _addchemicalpotential_lattice_bdg!(H, lat, μ; kwargs...)
Operators.addchemicalpotential!(H::BdGOperator, lat::Lattices.Lattice, μ::Function; kwargs...) =
    _addchemicalpotential_lattice_bdg!(H, lat, μ; kwargs...)
Operators.addchemicalpotential!(H::BdGOperator, lat::Lattices.Lattice, μ; kwargs...) =
    _addchemicalpotential_lattice_bdg!(H, lat, μ; kwargs...)

function electron(H::BdGOperator)
    d = hopdim(H)
    P = spzeros(ComplexF64, 2d, 2d)
    for i in 1:d
        P[i, i] = 1
    end
    Hops(Dict(zerokey(H) => P))
end



function Operators.localdensity(ρ::BdGOperator, lat::Lattices.Lattice)
    Operators.localdensity(Utils.getelectronsector(ρ), lat)
end

function Operators.localobservables(ρ::BdGOperator, lat::Lattices.Lattice)
    ρ0 = Utils.getelectronsector(ρ)
    M = Operators.localobservables(ρ0, lat)

    Δ= Superconductivity.getpairingsector(ρ) * (1im*Operators.getoperator(lat, "sy"))
    SC = Operators.localobservables(Δ, lat)

    M, SC
end

function Operators.expval(ρ::BdGOperator, args...; kwargs...)
    ρ0 = Utils.getelectronsector(ρ)
    Δ = Superconductivity.getpairingsector(ρ)
    lat = _find_lattice(args)
    if lat !== nothing
        Δ = Δ * (1im * Operators.getoperator(lat, "sy"))
    end

    M = Operators.expval(ρ0, args...; kwargs...)
    SC = Operators.expval(Δ, args...; kwargs...)

    M, SC
end

_find_lattice(args) = nothing
function _find_lattice(args::Tuple)
    for arg in args
        arg isa Lattices.Lattice && return arg
    end
    nothing
end

# ---------------------------------------------------------------------------
# Chemical-potential and filling that take the FULL BdG into account.
#
# The generic `Spectrum.chemicalpotential(H, ks, filling)` falls through to
# `Spectrum.chemicalpotential(Utils.getelectronsector(H), ks, filling)` via
# multiple dispatch, which strips Δ and fits μ to the normal-state bands.
# That's a common BdG approximation — accurate to O(Δ²/W) — but it breaks
# down for strong pairing or for non-conventional Δ that significantly
# rearranges the electron sector.
#
# These two methods do the right thing: bisect on μ such that the electron
# trace of the FULL BdG density matrix equals N · filling.
# ---------------------------------------------------------------------------

import Roots
import LatticeQM.Spectrum

"""
    Spectrum.filling(H::BdGOperator, ks; μ=0.0, T=0.01, kwargs...)

Electron filling fraction at chemical potential `μ` for the BdG Hamiltonian
`H`, computed from the full BdG density matrix (i.e. accounting for any
pairing block Δ). This is `Tr(ρ_el(k=0)) / N` averaged over `ks` weights.

`multimode` defaults to `:serial` here for two reasons:

1. The BdG `getdensitymatrix!` mutates `H` in place via
   `addchemicalpotential!(H, ±μ)` to apply the chemical potential. The
   shift is bracketed around the k-loop (added before, undone after), so
   threaded reads are safe — but the bracket is exception-fragile if a
   worker throws. Serial keeps the contract simple.

2. Bisection in `chemicalpotential(::BdGOperator)` calls this in a tight
   loop (Roots.find_zero); the per-call threading overhead exceeds the
   per-call eigen cost on small (e.g. 144-k) grids. Pass
   `multimode=:multithreaded` or `:distributed` for larger workloads.
"""
function Spectrum.filling(H::BdGOperator, ks, μ::Real=0.0;
                          T::Real=0.01, multimode::Symbol=:serial, kwargs...)
    d = hopdim(H)
    ρ = BdGOperator(zero(Utils.copyelectronsector(H)))
    Operators.getdensitymatrix!(ρ, H, ks, μ; T=T, multimode=multimode, kwargs...)
    real(tr(view(ρ.h[TightBinding.zerokey(ρ.h)], 1:d, 1:d))) / d
end

function Spectrum.filling(H::BdGOperator, μ::Real; nk=10, kwargs...)
    ks = Structure.regulargrid(nk=nk^2)
    Spectrum.filling(H, ks, μ; kwargs...)
end

"""
    Spectrum.chemicalpotential(H::BdGOperator, ks, target::Real; T=0.01, ...)

Chemical potential at which the BdG Hamiltonian `H` reaches the target
electron filling. Bisects on `μ` using `Spectrum.filling(H::BdGOperator)`
above (i.e. the *full* BdG, including Δ). The default
`Spectrum.chemicalpotential(::Any, ...)` strips Δ before fitting, which is
fine for weak coupling but wrong for strong pairing — this method overrides
that for `BdGOperator` so the SCF self-consistently solves μ even with a
non-zero gap.
"""
function Spectrum.chemicalpotential(H::BdGOperator, ks, target::Real;
                                     T::Real=0.01,
                                     μ_lo::Real=-100.0, μ_hi::Real=100.0,
                                     kwargs...)
    # Bracket from the normal-state spectrum so bisection converges fast on
    # well-conditioned models. We use `Spectrum.bandmatrix` on the electron
    # sector only as a guide for the bracket; the actual μ is fit against
    # the full BdG.
    h_el = Utils.getelectronsector(H)
    bands = Spectrum.bandmatrix(h_el, ks; multimode=:serial, hidebar=true,
                                progress_label="BdG chemical potential bracket")[1]
    em, eM = extrema(bands)
    lo = max(μ_lo, em - 1.0)
    hi = min(μ_hi, eM + 1.0)
    f(μ) = Spectrum.filling(H, ks, μ; T=T, kwargs...) - target
    Roots.find_zero(f, (lo, hi))
end
