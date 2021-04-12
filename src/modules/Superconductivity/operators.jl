import ..Operators
using ..Operators: addchemicalpotential!

using ..TightBinding: Hops, zerokey, hopdim
using SparseArrays#: spzeros

import ..Spectrum

function Operators.addchemicalpotential!(H::BdGOperator, μ::Real)
    d = hopdim(H)
    addhops!(H, BdGOperator(Hops(zerokey(H.h)=>complex(spzeros(d,d)+μ*I))))
end

function electron(H::BdGOperator)
    d = hopdim(H)

    BdGOperator(Hops(Dict(zerokey(H)=>zero(first(values(Spectrum.getelectronsector(H))))+1.0*I)))
end

function Operators.localobservables(ρ::BdGOperator, lat)
    ρ0 = Spectrum.getelectronsector(ρ)
    M = Operators.localobservables(ρ0, lat)

    Δ= Superconductivity.getpairingsector(ρ) * (1im*Operators.getoperator(lat, "sy"))
    SC = Operators.localobservables(Δ, lat)

    M, SC
end