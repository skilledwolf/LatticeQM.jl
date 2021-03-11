using Base
import ..TightBinding
using ..TightBinding: Hops, AnyHops, DenseHops, SharedDenseHops, hopdim, addhops!, addhops, zerokey


mutable struct BdGOperator{T}
    h::T
end

function BdGOperator(h0::T) where T<:AnyHops

    h1 = Hops()

    for R=keys(h0)
        h1[R] = [h0[R] zero(h0[R]); zero(h0[R]) -conj.(h0[R])]
    end

    BdGOperator{T}(h1)
end

function BdGOperator(h0::T, Δ0::T) where T<:AnyHops

    addpairing!(BdGOperator(h0), Δ0)
end

function addpairing!(H::BdGOperator{T}, Δ::AnyHops) where T<:AnyHops
    @assert hopdim(H) == hopdim(Δ) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    Δ1 = Hops()

    for R=keys(Δ)
        Δ1[R] = [zero(Δ[R]) Δ[R]'; Δ[R] zero(Δ[R])]
    end

    addhops!(H,Δ1)
end

function addelectronsector!(H::BdGOperator{T}, h::AnyHops) where T<:AnyHops
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    H2 = Hops()
    for R=keys(h)
        H2[R] = [h[R] zero(h[R]); zero(h[R]) zero(h[R])]
    end

    addhops!(H,H2)
end

function addholesector!(H::BdGOperator{T}, h::AnyHops) where T<:AnyHops
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    H2 = Hops()
    for R=keys(h)
        H2[R] = [zero(h[R]) zero(h[R]); zero(h[R]) zero(h[R]) h[R]]
    end

    addhops!(H,H2)
end

function TightBinding.getelectronsector(H::BdGOperator{T}) where T<:AnyHops
    hops = Hops()
    d = hopdim(H)

    for R=keys(H)
        hops[R] = H[R][1:d,1:d] #select electron part (no pairing, no holes)
    end

    hops
end

function getpairingsector(H::BdGOperator{T}) where T<:AnyHops
    Δ1 = Hops()
    d = hopdim(H)

    for R=keys(H)
        Δ1[R] = H[R][1:d,d+1:2*d] #select electron part (no pairing, no holes)
    end

    Δ1
end


# Indexing interface

Base.values(H::BdGOperator) = Base.values(H.h)
Base.keys(H::BdGOperator{<:AnyHops}) = Base.keys(H.h)
Base.getindex(H::BdGOperator,i) = Base.getindex(H.h,i)
Base.setindex!(H::BdGOperator,v,i) = Base.setindex!(H.h,v,i)
Base.firstindex(H::BdGOperator) = Base.firstindex(H.h)
Base.lastindex(H::BdGOperator) = Base.lastindex(H.h)


# Interface to tight binding module

TightBinding.hopdim(H::BdGOperator{T}) where T<:AnyHops = div(TightBinding.hopdim(H.h),2)

function TightBinding.DenseHops(H::BdGOperator{T}) where T<:AnyHops
    h = TightBinding.DenseHops(H.h)
    BdGOperator{typeof(h)}(h)
end

function TightBinding.SharedDenseHops(H::BdGOperator{T}) where T<:AnyHops
    h = TightBinding.SharedDenseHops(H.h)
    BdGOperator{typeof(h)}(h)
end

TightBinding.addhops!(H1::BdGOperator{T1}, H2::BdGOperator{T2}) where {T1<:AnyHops, T2<:AnyHops} = TightBinding.addhops!(H1, H2.h)
function TightBinding.addhops!(H1::BdGOperator{T}, h::AnyHops) where T<:AnyHops
    TightBinding.addhops!(H1.h, h)
    H1
end

TightBinding.addhops(H1::BdGOperator{T1}, H2::BdGOperator{T2}) where {T1<:AnyHops, T2<:AnyHops} = TightBinding.addhops(H1, H2.h)
function TightBinding.addhops(H1::BdGOperator{T}, h::AnyHops) where T<:AnyHops
    Hnew = deepcopy(H1)
    TightBinding.addhops!(Hnew.h, h)
    Hnew
end

TightBinding.zerokey(H::BdGOperator) = TightBinding.zerokey(H.h)