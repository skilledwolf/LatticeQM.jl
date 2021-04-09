using Base
import ..TightBinding
using ..TightBinding: Hops, AnyHops, DenseHops, SharedDenseHops, hopdim, addhops!, addhops, zerokey


mutable struct BdGOperator{T}
    h::T
end

(H::BdGOperator)(k) = H.h(k) # This will make the type callable


function BdGOperator(h0::T) where T<:Hops

    h1 = TightBinding.empty(h0) #Hops()

    for R=keys(h0)
        h1[R] = [h0[R] zero(h0[R]); zero(h0[R]) -conj.(h0[R])]
    end

    BdGOperator{T}(h1)
end

function BdGOperator(h0::T, Δ0::T) where T<:Hops

    addpairing!(BdGOperator(h0), Δ0)
end

function addpairing!(H::BdGOperator{T}, Δ::Hops) where T<:Hops
    @assert hopdim(H) == hopdim(Δ) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    Δ1 = Hops()

    for R=keys(Δ)
        Δ1[R] = [zero(Δ[R]) Δ[R]; Δ[-R]' zero(Δ[R])]
    end

    addhops!(H,Δ1)
end

function addelectronsector!(H::BdGOperator{T}, h::Hops) where T<:Hops
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    H2 = Hops()
    for R=keys(h)
        H2[R] = [h[R] zero(h[R]); zero(h[R]) zero(h[R])]
    end

    addhops!(H,H2)
end

function addholesector!(H::BdGOperator{T}, h::Hops) where T<:Hops
    @assert hopdim(H) == hopdim(h) "Mismatch between matrix sizes of BdG Operator and pairing matrices."

    H2 = Hops()
    for R=keys(h)
        H2[R] = [zero(h[R]) zero(h[R]); zero(h[R]) zero(h[R]) h[R]]
    end

    addhops!(H,H2)
end

# Interface to Spectrum
import ..Spectrum

function Spectrum.getelectronsector(H::BdGOperator{T}) where T<:Hops
    hops = Hops()
    d = hopdim(H)

    for R=keys(H)
        hops[R] = H[R][1:d,1:d] #select electron part (no pairing, no holes)
    end

    hops
end

function getpairingsector(H::BdGOperator{T}) where T<:Hops
    Δ1 = Hops()
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

function TightBinding.DenseHops(H::BdGOperator)
    h = TightBinding.DenseHops(H.h)
    BdGOperator{typeof(h)}(h)
end

function TightBinding.SharedDenseHops(H::BdGOperator)
    h = TightBinding.SharedDenseHops(H.h)
    BdGOperator{typeof(h)}(h)
end

TightBinding.addhops!(H1::BdGOperator, H2::BdGOperator) = TightBinding.addhops!(H1, H2.h)
function TightBinding.addhops!(H1::BdGOperator, h::Hops)
    TightBinding.addhops!(H1.h, h)
    H1
end

TightBinding.addhops(H1::BdGOperator, H2::BdGOperator) = TightBinding.addhops(H1, H2.h)
function TightBinding.addhops(H1::BdGOperator, h::Hops)
    Hnew = deepcopy(H1)
    TightBinding.addhops!(Hnew.h, h)
    Hnew
end

TightBinding.zerokey(H::BdGOperator) = TightBinding.zerokey(H.h)