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
    N = length(first(keys(hoppings)))
    1 <= d <= N || throw(ArgumentError("Fourier dimension must lie in 1:$N"))
    newhops = Dict{Vector{Int},AbstractMatrix{ComplexF64}}()
    for (key, h) in hoppings
        reduced = key[(1:N) .!= d]
        block = h * cis(2π * k * key[d])
        newhops[reduced] = haskey(newhops, reduced) ? newhops[reduced] + block : block
    end

    # Infer a concrete matrix type when all resulting blocks share one.
    Dict((R => block for (R, block) in newhops)...)
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

import SparseArrays
# Constrain the K type parameter explicitly so the SparseHops method is
# unambiguously more specific than the generic Hops method on Julia ≥ 1.12.
SparseArrays.issparse(hops::Hops{K,T}) where {K,T} = false
SparseArrays.issparse(hops::SparseHops{K}) where {K<:AbstractVector{<:Int}} = true

# (H::Hops{K,T})(k) where {K,T} = (h0::T = fouriersum(H, k); h0) #Hermitian(fouriersum(H, k))
(H::Hops{K,T})(k) where {K,T} = fouriersum(H.data, k)
fouriersum!(out::AbstractMatrix, H::Hops{K,T}, k) where {K,T} = fouriersum!(out, H.data, k)

fouriersum(hoppings::Hops, k::Real, d::Int) = Hops(fouriersum(hoppings.data, k, d))

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

# Zero-allocation Bloch matrix construction for any AbstractHops. This is the
# fast path that bandmatrix and downstream k-loops dispatch into; the generic
# Spectrum._build_H! fallback (which copies from a freshly-allocated H(k)) is
# only used for caller-supplied operator types that aren't AbstractHops.
Spectrum._build_H!(out::AbstractMatrix, H::AbstractHops, k) = fouriersum!(out, H, k)


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
    # Snapshot the keys: the `else` branch inserts H[-R] while iterating,
    # and mutating a Dict during iteration is undefined behavior (skipped or
    # repeated keys after a rehash) — exactly in the unpaired-keys case this
    # function exists to repair.
    Hkeys = collect(keys(H))
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
function addhops!(hops::Hops, newhops::Hops...)
    for other in newhops, (R, block) in other
        if haskey(hops, R)
            hops[R] = hops[R] + block
        else
            hops[copy(R)] = copy(block)
        end
    end
    hops
end
addhops(hops::Hops, newhops::Hops...) = (H=deepcopy(hops); addhops!(H,newhops...)) #merge(+, hops, newhops...)

Base.:*(h::Hops, s::Number) = multiplyhops(h,s)
Base.:*(s::Number, h::Hops) = multiplyhops(h,s)
Base.:*(h1::Hops, h2::Hops) = multiplyhops(h1,h2)
Base.:*(h1::Hops, h2::AbstractMatrix) = multiplyhops(h1,h2)
Base.:*(h1::AbstractMatrix, h2::Hops) = multiplyhops(h1,h2)
multiplyhops(M::AbstractMatrix, h::Hops) = Hops(copy(R) => M * block for (R, block) in h)
multiplyhops(h::Hops, M::AbstractMatrix) = Hops(copy(R) => block * M for (R, block) in h)
function multiplyhops(h1::Hops, h2::Hops)
    out = Hops()
    for (R, A) in h1, (S, B) in h2
        key = R + S
        block = A * B
        out[key] = haskey(out, key) ? out[key] + block : block
    end
    isempty(out.data) ? out : Hops(R => block for (R, block) in out)
end
multiplyhops(h::Hops, s::Number) = Hops(k=>h[k]*s for k=keys(h))
multiplyhops(s::Number, h::Hops) = multiplyhops(h, s)

"""
Combine the linear spaces of two hopping models without modifying either input.
"""
function directsum(h1::Hops, h2::Hops)
    d1 = hopdim(h1)
    d2 = hopdim(h2)
    D = d1+d2

    out = SparseHops()
    for δL in union(keys(h1), keys(h2))
        tmp = spzeros(ComplexF64, D, D)
        haskey(h1, δL) && (tmp[1:d1, 1:d1] = h1[δL])
        haskey(h2, δL) && (tmp[d1+1:D, d1+1:D] = h2[δL])
        out[copy(δL)] = tmp
    end

    out
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


function efficientzero(ρ::Hops)
    L = length(ρ)
    @assert L > 0 "Must have at least one hopping element."

    dims = size(first(values(ρ)))
    
    A = zeros(ComplexF64, dims..., L) #zeros(eltype(valtype(ρ)), dims..., L)

    A, collect(keys(ρ))
end

function flexibleformat!(ρ::Hops, A::AbstractArray, keylist::AbstractVector)
    for (j_,L)=enumerate(keylist)
        # ρ[L][:,:] .= m[:,:]
        # copyto!(ρ[L][:,:], A[:,:,j_])
        ρ[L] .= A[:,:,j_]
    end
    ρ
end
