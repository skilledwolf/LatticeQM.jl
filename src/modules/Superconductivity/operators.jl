import ..Operators
using ..Operators: addchemicalpotential!

using ..TightBinding: Hops, zerokey, hopdim
using SparseArrays#: spzeros

import ..Spectrum

Spectrum.dim(h::BdGOperator, args...) = Spectrum.dim(h.h, args...)

function Operators.addchemicalpotential!(H::BdGOperator, μ::Real)
    d = hopdim(H)
    addhops!(H, BdGOperator(Hops(zerokey(H.h)=>complex(spzeros(d,d)+μ*I))))
end

import ..Structure: Lattices

function Operators.addchemicalpotential!(H::BdGOperator, lat::Lattices.Lattice, μ; kwargs...)
    d = hopdim(H)

    chemhops = Hops(zerokey(H)=>zeros(ComplexF64,d,d))
    Operators.addchemicalpotential!(chemhops, lat, μ; kwargs...)

    addhops!(H, BdGOperator(chemhops))
end

function electron(H::BdGOperator)
    d = hopdim(H)
    BdGOperator(Hops(Dict(zerokey(H)=>zero(first(values(Spectrum.getelectronsector(H))))+1.0*I)))
end



function Operators.localdensity(ρ::BdGOperator, lat::Lattices.Lattice)
    Operators.localdensity(Spectrum.getelectronsector(ρ), lat)
end

function Operators.localobservables(ρ::BdGOperator, lat::Lattices.Lattice)
    ρ0 = Spectrum.getelectronsector(ρ)
    M = Operators.localobservables(ρ0, lat)

    Δ= Superconductivity.getpairingsector(ρ) * (1im*Operators.getoperator(lat, "sy"))
    SC = Operators.localobservables(Δ, lat)

    M, SC
end

function Operators.expval(ρ::BdGOperator, args...; kwargs...)
    ρ0 = Spectrum.getelectronsector(ρ)
    Δ= Superconductivity.getpairingsector(ρ) * (1im*Operators.getoperator(lat, "sy"))

    M = Operators.expval(ρ0, args...; kwargs...)
    SC = Operators.expval(Δ, args...; kwargs...)

    M, SC
end

