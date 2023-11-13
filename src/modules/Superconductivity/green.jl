
using ..Operators: densitymatrix, densitymatrix!
using ..TightBinding: zerokey


function Operators.getdensitymatrix!(ρ::BdGOperator{<:AnyHops}, H::BdGOperator{<:AnyHops}, ks::AbstractMatrix, μ::Real=0; kwargs...)
    d = hopdim(H)

    H[zerokey(H)][d+1:2*d,d+1:2*d] += 2*μ*I
    res = Operators.getdensitymatrix!(ρ.h, H.h, ks, μ; kwargs...)
    H[zerokey(H)][d+1:2*d,d+1:2*d] -= 2*μ*I

    res
end