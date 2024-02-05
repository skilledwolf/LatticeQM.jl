######################################################################################
# Bloch construction
######################################################################################

using LinearAlgebra: dot
using LinearAlgebra: Hermitian

fourierphase(k,δL) = exp(1.0im * 2π * dot(k, δL))

atleast1d(x::Number) = [x]
atleast1d(x::AbstractArray) = x

function fouriersum(hoppings, k)
    k=atleast1d(k)
    sum(t .* fourierphase(k, δL) for (δL, t) in hoppings)
end

function fouriersum!(out::AbstractMatrix, hoppings, k)
    out .= 0
    k=atleast1d(k)
    for (δL, t) in hoppings
        out .+= t .* fourierphase(k, δL)
    end
    out
end

function fouriersum(hoppings::Dict{K,T}, k::Real, d::Int) where {K,T}
    N=length(zerokey(hoppings))
    newhops = Dict()

    for key in keys(hoppings)
        newhops[key[(1:N).!=d]] = zeros(typeof(hoppings[key]), size(hoppings[key]))
    end

    for (key,h)=hoppings
        newhops[key[(1:N).!=d]] .+= h .* fourierphase([k],[key[d]])
    end

    newhops
end

######################################################################################
# Tight binding operator definitions
######################################################################################

import SharedArrays: SharedArray, SharedMatrix
import SparseArrays
import SparseArrays: SparseMatrixCSC, sparse, spzeros
import LatticeQM.Utils
# import ..Utils: dense

abstract type AbstractHops{K,T} end

struct Hops{K<:AbstractVector{<:Int}, T<:AbstractMatrix{ComplexF64}} <: AbstractHops{K,T}
    data::Dict{K,T}
end

const DenseHops{K}       = Hops{K,Matrix{ComplexF64}}
const SharedDenseHops{K} = Hops{K,SharedMatrix{ComplexF64}}
const SparseHops{K}      = Hops{K,SparseMatrixCSC{ComplexF64,Int64}}
const SubarrayHops{K} = Hops{K,<:SubArray} # todo: run checks

gethopsview(h::Hops) = Hops(Dict(L => view(M, :, :) for (L, M) in h))
gethopsview(h::SubarrayHops) = h

import SparseArrays
SparseArrays.issparse(hops::Hops) = false
SparseArrays.issparse(hops::SparseHops) = true 

# (H::Hops{K,T})(k) where {K,T} = (h0::T = fouriersum(H, k); h0) #Hermitian(fouriersum(H, k))
(H::Hops{K,T})(k) where {K,T} = fouriersum(H.data, k)
fouriersum!(out::AbstractMatrix, H::Hops{K,T}, k) where {K,T} = fouriersum!(out, H.data, k)

fouriersum(hoppings::Hops{K,T}, k::Real, d::Int) where {K,T} = Hops{K,T}(fouriersum(hoppings.data, k, d))

# Hops{K,T}(d::AbstractDict) where {K<:AbstractVector{Int},T<:AbstractMatrix{ComplexF64}} = (print("miau"); ) # Hops{K,T}(Dict{K,T}(d...))

# todo: change the logic here:
# Hops is very flexible, but DenseHops, SparseHops, and SharedDenseHops create unnecessary copies of Hops first

Hops(hops::Hops) = hops
Hops(T::Type=AbstractMatrix{ComplexF64}, K::Type=Vector{Int}) = Hops{K,T}(Dict{K,T}()) #
Hops(M::AbstractMatrix, d::Int=2) = Hops(Dict(zeros(Int,d)=>M))
Hops(kv::Pair...) = Hops(Dict(k=>v for (k,v) in kv))
Hops(G::Base.Generator) = Hops(Dict(G...))

# Dense conversion
Utils.dense(hops::DenseHops) = hops # does not create new instance, pass by reference
Utils.dense(hops::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}} = DenseHops{K}(Dict(k => Utils.dense(v) for (k, v) in hops.data)) # creates new instances, copies data
Utils.densecopy(hops::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}} = DenseHops{K}(Dict(k => Utils.densecopy(v) for (k, v) in hops.data))
DenseHops(args...; kwargs...) = Utils.dense(Hops(args...; kwargs...)) # uses constructors of Hops to create object, and dense(...) to potentially create yet another copy. This seems wasteful

# Sparse conversion
SparseArrays.sparse(hops::SparseHops) = hops
SparseArrays.sparse(hops::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}} = SparseHops{K}(Dict(k => sparse(v) for (k, v) in hops.data))
SparseHops(args...; kwargs...) = sparse(Hops(args...; kwargs...)) # uses constructors of Hops to create object, and sparse(...) to potentially create yet another copy. This seems wasteful

# Shared conversion
shareddense(hops::SharedDenseHops) = hops
shareddense(hops::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}} = SharedDenseHops{K}(Dict(k => SharedArray(Matrix{ComplexF64}(v)) for (k, v) in hops.data))
SharedDenseHops(args...; kwargs...) = shareddense(Hops(args...; kwargs...))

import LatticeQM.Utils
autoconversion(hops::Hops, N::Int) = (N < MAX_DENSE + 1) ? Utils.dense(hops) : sparse(hops)
autoconversion(hops::Hops, N::Int, type::Symbol) = (type == :auto) ? autoconversion(hops, N) : (type == :sparse) ? sparse(hops) : Utils.dense(hops)

Utils.getelectronsector(H::Hops) = H
Utils.copyelectronsector(H::Hops) = deepcopy(H)

import LatticeQM.Spectrum 
Spectrum.sanatize_distributed_hamiltonian(H::DenseHops) = shareddense(H)

# Size
Base.size(H::Hops, args...) = Base.size(first(values(H.data)), args...)
Base.broadcastable(H::Hops) = Ref(H)

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
Base.empty(H::Hops{K,T}) where {K,T} = Hops{K,T}(Dict{K,T}())

function Base.fill!(H::Hops, x)
    for δL in keys(H)
        fill!(H[δL], x)
    end
    H
end

zero_matrix(::Type{Hops{K,T}}, N1, N2) where {K,T} = zeros(ComplexF64, N1, N2)
zero_matrix(::Type{SparseHops{K}}, N1, N2) where {K} = spzeros(ComplexF64, N1, N2)
zero_matrix(::Type{SharedDenseHops{K}}, N1, N2) where {K} = SharedMatrix(zeros(ComplexF64, N1, N2))

function Base.zero(h::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}}
    Hops{K,T}(Dict{K,T}(δL => zero(h[δL]) for δL in keys(h)))
end

function zero_dense(h::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}}
    d = Dict{K,Matrix{ComplexF64}}(δL => zeros(ComplexF64, size(h[δL])) for δL in keys(h))
    DenseHops{K}(d)
end

function zero_sparse(h::Hops{K,T}) where {K,T<:AbstractMatrix{ComplexF64}}
    d = Dict{K,SparseMatrixCSC{ComplexF64,Int64}}(δL => spzeros(ComplexF64, size(h[δL])) for δL in keys(h))
    SparseHops{K}(d)
end

function ishermitian(H::Hops; tol=sqrt(eps()))
    for R=keys(H)
        if !(haskey(H,-R) && maximum(abs.(H[R].-H[-R]'))<tol) #"H[R] does not have H[-R] partner."
            return false   
        end
    end
    return true
end

function hermitianize!(H::Hops) # hermitian means H[R] = H[-R]'
    Hkeys = keys(H)
    for R in Hkeys
        if haskey(H, -R)
            @. H[R] = (H[R] + H[-R]')/2
            @. H[-R] = H[R]'
        else 
            H[-R] = H[R]'
        end
    end
    H
end

zerokey(h::Hops) = zero(first(keys(h)))
getzero(h::Hops) = h[zerokey(h)]
setzero!(h::Hops, M::AbstractMatrix) = (h[zerokey(h)].=M; h)

hopdim(hops::Hops) = size(hops,1)

Base.:+(h1::Hops, h2::Hops) = addhops(h1,h2)
Base.:-(h1::Hops, h2::Hops) = addhops(h1,(-1)*h2)
addhops!(hops::Hops, newhops::Hops...) = (mergewith!(+, hops.data, map(x->x.data,newhops)...); hops)
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

const MAX_DENSE = 500::Int
const MAX_DIAGS = 100::Int

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

function trim!(ρ::Hops{K,<:AbstractMatrix{ComplexF64}}; kwargs...) where {K}
    for δL in keys(ρ)
        if all(isapprox.(ρ[δL], 0; kwargs...))
            delete!(ρ.data, δL)
        end
    end
    ρ
end
function trim!(ρ::SparseHops{K}; kwargs...) where {K}
    for δL in keys(ρ)
        ρ[δL] = SparseArrays.dropzeros!(ρ[δL])
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
        A[:,:,i] .= ρ[δL]
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
    Dict(L=>m for (L,m)=zip(keylist, eachslice(A; dims=3)))
end

function flexibleformat!(ρ::Hops, A::AbstractArray, keylist::AbstractVector)
    for (j_,L)=enumerate(keylist)
        # ρ[L][:,:] .= m[:,:]
        # copyto!(ρ[L][:,:], A[:,:,j_])
        ρ[L] .= A[:,:,j_]
    end
    ρ
end