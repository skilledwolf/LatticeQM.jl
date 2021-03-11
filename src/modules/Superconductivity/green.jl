
# Interface to Green module 
import ..Green

using ..Green: densitymatrix, densitymatrix!
using ..TightBinding: zerokey

function Green.densitymatrix(H::BdGOperator{<:AnyHops}, args...; kwargs...)
    Green.densitymatrix(H.h, args...; kwargs...)
end

function Green.densitymatrix!(ρ::BdGOperator{<:AnyHops}, H::BdGOperator{<:AnyHops}, ks::AbstractMatrix, μ::Real=0; kwargs...)
    d = hopdim(H)

    H[zerokey(H)][d+1:2*d,d+1:2*d] += 2*μ*I
    res = Green.densitymatrix!(ρ.h, H.h, ks, μ; kwargs...)
    H[zerokey(H)][d+1:2*d,d+1:2*d] -= 2*μ*I

    res
end