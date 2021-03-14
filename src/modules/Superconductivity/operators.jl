import ..Operators
using ..Operators: addchemicalpotential!

using ..TightBinding: Hops, zerokey, hopdim
using SparseArrays: spzeros

function Operators.addchemicalpotential!(H::BdGOperator{<:AnyHops}, μ::Real)
    d = hopdim(H)
    addhops!(H, BdGOperator(Hops(zerokey(H.h)=>spzeros(d,d)+μ*I)))
end

function electron(H::BdGOperator{<:AnyHops})
    d = hopdim(H)

    BdGOperator(Hops(zerokey(H)=>zero(first(values(TightBinding.getelectronsector(H))))+1.0*I))
end