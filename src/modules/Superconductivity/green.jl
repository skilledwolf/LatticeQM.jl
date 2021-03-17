
# Interface to Green module 
import ..Green

using ..Green: densitymatrix, densitymatrix!
using ..TightBinding: zerokey

function Green.getdensitymatrix(H::BdGOperator{<:AnyHops}, args...; kwargs...)
    Green.getdensitymatrix(H.h, args...; kwargs...)
end

function Green.getdensitymatrix!(ρ::BdGOperator{<:AnyHops}, H::BdGOperator{<:AnyHops}, ks::AbstractMatrix, μ::Real=0; kwargs...)
    d = hopdim(H)

    H[zerokey(H)][d+1:2*d,d+1:2*d] += 2*μ*I
    res = Green.getdensitymatrix!(ρ.h, H.h, ks, μ; kwargs...)
    H[zerokey(H)][d+1:2*d,d+1:2*d] -= 2*μ*I

    res
end