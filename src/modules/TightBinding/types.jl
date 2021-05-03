######################################################################################
# Bloch construction
######################################################################################

using LinearAlgebra: dot

fourierphase(k,δL) = exp(1.0im * 2π * dot(k, δL))

fouriersum(hoppings, k::AbstractVector) = sum(t .* fourierphase(k, δL) for (δL,t) in hoppings)
fouriersum(hoppings, k::Float64) = sum(t .* fourierphase([k], δL) for (δL,t) in hoppings)

function fouriersum(hoppings, k::Real, d::Int)
    N=length(zerokey(hoppings))

    K = unique(collect(key[(1:N).!=d] for key=keys(hoppings)))
    newhops = Dict(key => zero(getzero(hoppings)) for key in K)

    for (key,h)=hoppings
        newhops[key[(1:N).!=d]] .+= h .* fourierphase([k],[key[d]])
    end

    newhops
end

function getbloch(hoppings,args...)
    function h(k)
        fouriersum(hoppings, k,args...)
    end
    h
end

######################################################################################
# Tight binding operator definitions
######################################################################################

import SharedArrays: SharedArray
import SparseArrays: SparseMatrixCSC, sparse, spzeros

mutable struct Hops{T<:AbstractMatrix{ComplexF64}}
    data::Dict{Vector{Int}, T}
end

const AnyHops = Hops#{AbstractMatrix{ComplexF64}}

(H::Hops)(k) = fouriersum(H,k) # This will make the type callable

Hops(T=AbstractMatrix{ComplexF64}) = Hops(Dict{Vector{Int},T}())
Hops(M::AbstractMatrix, d::Int=2) = Hops(Dict(zeros(Int,d)=>M))
Hops(kv::Pair...) = Hops(Dict(k=>v for (k,v) in kv))
Hops(d::AbstractDict) = Hops(Dict(d...))
Hops(G::Base.Generator) = Hops(Dict(G...))

DenseHops() = Hops(Matrix{ComplexF64})
DenseHops(H::Hops) = DenseHops(H.data)
DenseHops(kv::Pair...) = Hops{Matrix{ComplexF64}}(Dict(k=>complex(Matrix(v)) for (k,v) in kv))
DenseHops(d::AbstractDict) = DenseHops(d...)
DenseHops(G::Base.Generator) = DenseHops(Dict(G...))

SharedDenseHops(H::Hops) = SharedDenseHops(H.data)
SharedDenseHops(kv::Pair...) = Hops{Matrix{ComplexF64}}(Dict(k=>SharedArray(Matrix(complex(v))) for (k,v) in kv))
SharedDenseHops(d::AbstractDict) = SharedDenseHops(d...)
SharedDenseHops(G::Base.Generator) = SharedDenseHops(Dict(G...))

SparseHops() = Hops(SparseMatrixCSC{Complex{Float64},Int64})
SparseHops(H::Hops) = SparseHops(H.data)
SparseHops(kv::Pair...) = Hops{SparseMatrixCSC{ComplexF64,Int64}}(Dict(k=>sparse(complex(v)) for (k,v) in kv))
SparseHops(d::AbstractDict) = SparseHops(d...)
SparseHops(G::Base.Generator) = SparseHops(Dict(G...))

# Interface to spectrum
import ..Spectrum
Spectrum.getelectronsector(H::Hops) = H

# Size
Base.size(H::Hops, args...) = Base.size(first(values(H.data)), args...)

# Interface for iteration and item access
Base.get(H::Hops, args...) = Base.get(H.data, args...)
Base.haskey(H::Hops, args...) = Base.haskey(H.data, args...)
Base.iterate(H::Hops) = Base.iterate(H.data)
Base.iterate(H::Hops, state) = Base.iterate(H.data, state)
Base.length(H::Hops) = Base.length(H.data)
Base.eltype(H::Hops) = Base.eltype(H.data)

Base.values(H::Hops) = Base.values(H.data)
Base.keys(H::Hops) = Base.keys(H.data)
Base.firstindex(H::Hops) = Base.firstindex(H.data)
Base.lastindex(H::Hops) = Base.lastindex(H.data)
Base.getindex(H::Hops,i) = Base.getindex(H.data,i)
Base.setindex!(H::Hops,v,i) = Base.setindex!(H.data,v,i)

Base.empty!(H::Hops) = (Base.empty!(H.data); H)
Base.empty(H::Hops) = (H2=Base.empty!(deepcopy(H)); H2)


function Base.zero(h::Hops; format=:auto)
    if format==:dense
        return zero_dense(h)
    elseif format==:sparse
        return zero_sparse(h)
    else
        return zero_auto(h)
    end
end
function zero_dense(h::Hops)
    ρ = Hops()
    for δL=keys(h)
        ρ[δL] = zeros(ComplexF64, size(h[δL]))
    end
    ρ
end
function zero_sparse(h::Hops)
    ρ = Hops()
    for δL=keys(h)
        ρ[δL] = spzeros(ComplexF64, size(h[δL]))
    end
    ρ
end
function zero_auto(h::Hops)
    ρ = Hops()
    for δL=keys(h)
        ρ[δL] = zero(h[δL])
    end
    ρ
end


function ishermitian(H::Hops; tol=sqrt(eps()))
    for R=keys(H)
        if !(haskey(H,-R) && maximum(abs.(H[R].-H[-R]'))<tol) #"H[R] does not have H[-R] partner."
            return false   
        end
    end

    return true
end

getelectronsector(H::Function) = H
getelectronsector(H::AbstractMatrix) = H
getelectronsector(H::Hops) = H

zerokey(h::Hops) = zero(first(keys(h)))
getzero(h::Hops) = h[zerokey(h)]
setzero!(h::Hops, M::AbstractMatrix) = (h[zerokey(h)].=M; h)

hopdim(hops::Hops) = size(hops,1)

Base.:+(h1::Hops, h2::Hops) = addhops(h1,h2)
Base.:-(h1::Hops, h2::Hops) = addhops(h1,(-1)*h2)
addhops!(hops::Hops, newhops::Hops...) = (merge!(+, hops.data, map(x->x.data,newhops)...); hops)
addhops(hops::Hops, newhops::Hops...) = (H=deepcopy(hops); addhops!(H,newhops...)) #merge(+, hops, newhops...)

Base.:*(h::Hops, s::Number) = multiplyhops(h,s)
Base.:*(s::Number, h::Hops) = multiplyhops(h,s)
Base.:*(h1::Hops, h2::Hops) = multiplyhops(h1,h2)
Base.:*(h1::Hops, h2::AbstractMatrix) = multiplyhops(h1,h2)
Base.:*(h1::AbstractMatrix, h2::Hops) = multiplyhops(h1,h2)
multiplyhops(h1::AbstractMatrix, h2::Hops) = multiplyhops(Hops(h1),h2)
multiplyhops(h1::Hops, h2::AbstractMatrix) = multiplyhops(h1,Hops(h2))
multiplyhops(h1::Hops, h2::Hops) = Hops(k=>h1[k]*h2[k] for k=intersect(keys(h1),keys(h2)))
multiplyhops(h::Hops, s::Number) = Hops(k=>h[k]*s for k=keys(h))
multiplyhops(s::Number, h::Hops) = multiplyhops(h, s)

"""
Naive implementation of combining the linear spaces of two hopping models.
"""
function directsum(h1::Hops, h2::Hops)
    d1 = hopdim(h1)
    d2 = hopdim(h2)
    D = d1+d2

    for δL in keys(h1)
        tmp = spzeros(D,D)
        tmp[1:d1,1:d1] .= h1[δL]
        h1[δL] = tmp
    end
    for δL in keys(h2)
        tmp = spzeros(D,D)
        tmp[d1+1:D,d1+1:D] .= h2[δL]
        h2[δL] = tmp
    end

    addhops(h1,h2)
end

function Base.kron(a, b::Hops)
    b2 = deepcopy(b)

    for (δL, t) in b
        b2[δL] = kron(a, t) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    b2
end

function Base.kron(a::Hops, b)
    a2 = deepcopy(a)

    for (δL, t) in a
        a2[δL] = kron(t, b) # add spin degree of freedom # we made the choice to group the matrix in spin blocks
    end

    a2
end


import ..Utils

function addspin(hoppings, mode=:nospin) #::AbstractHops
    if mode==:nospin || mode==:id
        return hoppings
    elseif mode==:spinhalf || mode==:σ0
        hoppings = kron(hoppings, Utils.σ0)
    elseif mode==:σx
        hoppings = kron(hoppings, Utils.σX)
    else
        error("Do not recognize mode '$mode' in addspin(...).")
    end

    hoppings
end

const MAX_DENSE = 500
const MAX_DIAGS = 100

decidetype(hops::Hops) = decidetype(hopdim(hops))

function decidetype(N::Int)
    if N < MAX_DENSE + 1
        return :dense
    else
        return :sparse
    end
end


function ensuretype(hops::Hops, format=:auto)
    if format==:auto
        # format = decidetype(hops) # old behaviour
        return hops
    end 

    if format==:dense
        hops = DenseHops(hops)
    elseif format==:sparse
        hops = SparseHops(hops)
    end

    hops
end

import LinearAlgebra: norm
import SparseArrays

function trim!(ρ::Hops; kwargs...)
    for δL in keys(ρ)
        if SparseArrays.issparse(ρ[δL])
            ρ[δL] = sparse(ρ[δL])
            SparseArrays.dropzeros!(ρ[δL])
        end
        if all(isapprox.(ρ[δL], 0; kwargs...))
            delete!(ρ.data, δL)
        end
    end
    ρ
end

function efficientformat(ρ::Hops)
    L = length(ρ)
    @assert L > 0 "Must have at least one hopping element."

    dims = size(first(values(ρ)))
    
    A = Array{ComplexF64}(undef, dims..., L) #Array{eltype(valtype(ρ))}(undef, dims..., L)
    
    keylist = []
    for (i,δL) in enumerate(keys(ρ))
        A[:,:,i] .= ρ[δL][:,:]
        append!(keylist, [δL])
    end
    
    A, keylist
end

function efficientzero(ρ::Hops)
    L = length(ρ)
    @assert L > 0 "Must have at least one hopping element."

    dims = size(first(values(ρ)))
    
    A = zeros(ComplexF64, dims..., L) #zeros(eltype(valtype(ρ)), dims..., L)

    A, collect(keys(ρ))
end

function flexibleformat(A::AbstractArray, keylist::AbstractVector)
    Dict(L=>Matrix(m) for (L,m)=zip(keylist, eachslice(A; dims=3)))
end

function flexibleformat!(ρ::Hops, A::AbstractArray, keylist::AbstractVector)
    for (j_,L)=enumerate(keylist)
        # ρ[L][:,:] .= m[:,:]
        # copyto!(ρ[L][:,:], A[:,:,j_])
        ρ[L][:,:] .= A[:,:,j_]
    end
    ρ
end