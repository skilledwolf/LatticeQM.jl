import LatticeQM.TightBinding
using LatticeQM.TightBinding: Hops, hopdim, addhops!, addhops, zerokey


mutable struct BdGOperator{T}
    h::T
end

(H::BdGOperator)(k) = H.h(k) # This will make the type callable, because T is callable


function BdGOperator(h0::T) where T<:Hops
    h1 = TightBinding.empty(h0) # create empty copy
    for R=keys(h0)
        h1[R] = [h0[R]         zero(h0[R]); 
                 zero(h0[R]) -conj.(h0[R])]
    end
    BdGOperator{T}(h1)
end

function BdGOperator(h0::T, Δ0::T) where {T<:Hops}
    H = BdGOperator(h0)
    addpairing!(H, Δ0)
    H
end

function addpairing!(H::BdGOperator{T}, Δ::T2) where {T<:Hops, T2<:Hops}
    @assert hopdim(H) == hopdim(Δ) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    Δ1 = T2(Dict())
    for R=keys(Δ)
        Δ1[R] = [zero(Δ[R])   Δ[R]; 
                 Δ[-R]'       zero(Δ[R])]
    end

    addhops!(H,Δ1) # effectively mergewith!(+,)
end

function addelectronsector!(H::BdGOperator{T}, h::T2) where {T<:Hops,T2<:Hops}
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    H2 = T2(Dict())
    for R=keys(h)
        H2[R] = [h[R] zero(h[R]); zero(h[R]) -h[R]]
    end
    # H2 = Hops(Dict(R=>[h[R] zero(h[R]); zero(h[R]) -h[R]] for R=keys(h)))
    addhops!(H,H2)
end

function getelectronview(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)

    # hops = Hops()
    # for R = keys(H)
    #     hops[R] = view(H[R], 1:d, 1:d) #select electron part (no pairing, no holes)
    # end
    # hops
    Hops(Dict(R => view(H[R], 1:d, 1:d) for R = keys(H)))
end

function getpairingview(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)

    # Δ1 = Hops()
    # for R = keys(H)
    #     Δ1[R] = view(H[R], 1:d, d+1:2*d) #select electron part (no pairing, no holes)
    # end
    # Δ1
    Hops(Dict(R => view(H[R], 1:d, d+1:2*d) for R = keys(H)))
end


import LatticeQM.Utils

function Utils.getelectronsector(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)
    h = Dict(R=>view(H[R], 1:d, 1:d) for R=keys(H))
    Hops(h) # creates a view into the original object # TODO: check if this is works as intended
end
function Utils.copyelectronsector(H::BdGOperator{T}) where {T<:Hops}
    d = hopdim(H)
    Hops(Dict(R => H[R][1:d, 1:d] for R = keys(H))) # creates a copy of the original object # TODO: check if this is works as intended
end

function getpairingsector(H::BdGOperator{T}) where T<:Hops
    Δ1 = T(Dict())
    d = hopdim(H)

    for R=keys(H)
        Δ1[R] = H[R][1:d, d+1:2*d] #select electron part (no pairing, no holes)
    end

    Δ1
end


# Indexing interface
Base.values(H::BdGOperator) = Base.values(H.h)
Base.keys(H::BdGOperator) = Base.keys(H.h)
Base.getindex(H::BdGOperator,i) = Base.getindex(H.h,i)
Base.setindex!(H::BdGOperator,v,i) = Base.setindex!(H.h,v,i)
Base.firstindex(H::BdGOperator) = Base.firstindex(H.h)
Base.lastindex(H::BdGOperator) = Base.lastindex(H.h)


# Interface to tight binding module
TightBinding.hopdim(H::BdGOperator) = div(TightBinding.hopdim(H.h),2)

TightBinding.hermitianize!(H::BdGOperator) = TightBinding.hermitianize!(H.h)
TightBinding.ishermitian(H::BdGOperator) = TightBinding.ishermitian(H.h)

import ..Utils
import SparseArrays 
import ..TightBinding
Utils.dense(H::BdGOperator) = (h0 = Utils.dense(H.h); BdGOperator{typeof(h0)}(h0))
SparseArrays.sparse(H::BdGOperator) = (h0 = SparseArrays.sparse(H.h); BdGOperator{typeof(h0)}(h0))
TightBinding.shareddense(H::BdGOperator) = (h0 = TightBinding.shareddense(H.h); BdGOperator{typeof(h0)}(h0))

Base.:-(h1::BdGOperator, h2::BdGOperator) = h1.h - h2.h
Base.:-(h1::BdGOperator, h2::Hops) = h1.h - h2
Base.:+(h1::BdGOperator, h2::BdGOperator) = h1.h + h2.h
Base.:+(h1::BdGOperator, h2::Hops) = h1.h + h2
TightBinding.addhops!(H1::BdGOperator, H2::BdGOperator) = TightBinding.addhops!(H1, H2.h)
TightBinding.addhops!(H1::BdGOperator, h::Hops) = (TightBinding.addhops!(H1.h, h); H1)
    

TightBinding.addhops(H1::BdGOperator, H2::BdGOperator) = TightBinding.addhops(H1, H2.h)
function TightBinding.addhops(H1::BdGOperator, h::Hops)
    Hnew = deepcopy(H1)
    TightBinding.addhops!(Hnew.h, h)
    Hnew
end

TightBinding.zerokey(H::BdGOperator) = TightBinding.zerokey(H.h)