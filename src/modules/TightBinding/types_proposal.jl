const MAX_DENSE = 500::Int
const MAX_DIAGS = 100::Int

######################################################################################
# Bloch construction
######################################################################################

using LinearAlgebra: dot

fourierphase(k, δL) = exp(1.0im * 2π * dot(k, δL))

atleast1d(x::Number) = [x]
atleast1d(x::AbstractArray) = x
fouriersum(hoppings::Hops, k) = (k = atleast1d(k); sum(hoppings.hoppingmatrices[:, :, i] .* fourierphase(k, hoppings.basis.data[:, i]) for i in 1:size(hoppings.basis.data, 2)))

function fouriersum(hoppings::Hops, k::Real, d::Int)
    N = size(hoppings.basis.data, 1)
    newmatrices = copy(hoppings.hoppingmatrices)
    for i in 1:size(hoppings.basis.data, 2)
        if hoppings.basis.data[d, i] != 0
            newmatrices[:, :, i] .*= fourierphase([k], [hoppings.basis.data[d, i]])
        end
    end
    Hops(hoppings.basis, newmatrices)
end

######################################################################################
# Tight binding operator definitions
######################################################################################

import SharedArrays: SharedArray, SharedMatrix
import SparseArrays
import SparseArrays: SparseMatrixCSC, sparse, spzeros
import ..Utils
import ..Utils: dense

abstract type AbstractHops{K,T} end

abstract type AbstractHoppingBasis{T} end

struct HoppingBasis{T} <: AbstractHoppingBasis{T}
    data::Matrix{T} # columns are basis vectors, nr of rows is the dimension of the lattice
end

const AbstractHoppingMatrices = AbstractArray{ComplexF64,3}

struct Hops{K<:AbstractHoppingBasis,T<:AbstractHoppingMatrices} <: AbstractHops{K,T}
    basis::K
    hoppingmatrices::T # [:,:,i_] is the hopping matrix for the i-th basis vector
end

const DenseHops{K} = Hops{K,Matrix{ComplexF64}}
const SharedDenseHops{K} = Hops{K,SharedMatrix{ComplexF64}}
const SparseHops{K} = Hops{K,SparseMatrixCSC{ComplexF64,Int64}}

import SparseArrays
SparseArrays.issparse(hops::Hops) = false
SparseArrays.issparse(hops::SparseHops) = true

(H::Hops{K,T})(k) where {K,T} = fouriersum(H, k)

Hops(hops::Hops) = hops
Hops(basis::AbstractHoppingBasis, matrices::AbstractHoppingMatrices) = Hops(basis, matrices)
Hops(basis::Matrix, matrices::AbstractHoppingMatrices) = Hops(HoppingBasis(basis), matrices)

# Dense conversion
Utils.dense(hops::DenseHops) = hops
Utils.dense(hops::Hops{K,T}) where {K,T<:AbstractHoppingMatrices} = DenseHops{K}(HoppingBasis(hops.basis.data), convert(Matrix{ComplexF64}, hops.hoppingmatrices))

DenseHops(args...; kwargs...) = dense(Hops(args...; kwargs...))

# Sparse conversion
SparseArrays.sparse(hops::SparseHops) = hops
SparseArrays.sparse(hops::Hops{K,T}) where {K,T<:AbstractHoppingMatrices} = SparseHops{K}(HoppingBasis(hops.basis.data), SparseArrays.sparse(hops.hoppingmatrices))

SparseHops(args...; kwargs...) = sparse(Hops(args...; kwargs...))

# Shared conversion
shareddense(hops::SharedDenseHops) = hops
shareddense(hops::Hops{K,T}) where {K,T<:AbstractHoppingMatrices} = SharedDenseHops{K}(HoppingBasis(hops.basis.data), SharedArrays.SharedArray(convert(Matrix{ComplexF64}, hops.hoppingmatrices)))

SharedDenseHops(args...; kwargs...) = shareddense(Hops(args...; kwargs...))

# Autoconversion
autoconversion(hops::Hops, N::Int) = (N < MAX_DENSE + 1) ? dense(hops) : sparse(hops)
autoconversion(hops::Hops, N::Int, type::Symbol) = (type == :auto) ? autoconversion(hops, N) : (type == :sparse) ? sparse(hops) : dense(hops)

import LatticeQM.Utils
Utils.getelectronsector(H::Hops) = H

import LatticeQM.Spectrum
Spectrum.sanatize_distributed_hamiltonian(H::DenseHops) = shareddense(H)


# Size
Base.size(H::Hops{K,T}, args...) where {K,T} = Base.size(H.hoppingmatrices, args...)
Base.broadcastable(H::Hops) = Ref(H)

# Interface for iteration and item access
# Note: These methods may need further adaptation based on specific iteration requirements.
Base.iterate(H::Hops{K,T}) where {K,T} = zip(eachcol(H.basis.data), eachslice(H.hoppingmatrices, dims=3))
Base.length(H::Hops{K,T}) where {K,T} = size(H.basis.data, 2)
Base.eltype(H::Hops{K,T}) where {K,T} = Tuple{typeof(H.basis.data[:, 1]),eltype(H.hoppingmatrices)}

# Access individual hopping matrices
Base.getindex(H::Hops{K,T}, i::Int) where {K,T} = H.hoppingmatrices[:, :, i]
Base.setindex!(H::Hops{K,T}, v, i::Int) where {K,T} = (H.hoppingmatrices[:, :, i] = v; H)

# Operations on the whole set of hopping matrices
Base.empty!(H::Hops{K,T}) where {K,T} = (fill!(H.hoppingmatrices, 0); H)
Base.empty(H::Hops{K,T}) where {K,T} = Hops(K(), T(size(H.hoppingmatrices)...))

function Base.fill!(H::Hops{K,T}, x) where {K,T}
    fill!(H.hoppingmatrices, x)
    H
end

# Zero matrix constructors
zero_matrix(::Type{Hops{K,T}}, N1, N2) where {K,T} = zeros(ComplexF64, N1, N2)
zero_matrix(::Type{SparseHops{K}}, N1, N2) where {K} = spzeros(ComplexF64, N1, N2)
zero_matrix(::Type{SharedDenseHops{K}}, N1, N2) where {K} = SharedMatrix(zeros(ComplexF64, N1, N2))

# Zero Hops
function Base.zero(h::Hops{K,T}) where {K,T<:AbstractHoppingMatrices}
    zero_matrices = map(_ -> zero_matrix(T, size(h.hoppingmatrices, 1), size(h.hoppingmatrices, 2)), 1:size(h.basis.data, 2))
    Hops(h.basis, cat(zero_matrices..., dims=3))
end

# Zero dense Hops
function zero_dense(h::Hops{K,T}) where {K,T<:AbstractHoppingMatrices}
    DenseHops(K(), convert(Matrix{ComplexF64}, zero(h)))
end

# Zero sparse Hops
function zero_sparse(h::Hops{K,T}) where {K,T<:AbstractHoppingMatrices}
    SparseHops(K(), SparseArrays.sparse(zero(h)))
end

# Check if Hops is Hermitian
function ishermitian(H::Hops; tol=sqrt(eps()))
    for i = 1:size(H.basis.data, 2)
        if !(maximum(abs.(H.hoppingmatrices[:, :, i] - H.hoppingmatrices[:, :, i]')) < tol)
            return false
        end
    end
    return true
end

# Utility functions related to zero hopping
function zerokey(H::Hops)
    # Find the column in the basis data that contains only zeros
    findfirst(all(isequal(0), col) for col in eachcol(H.basis.data))
end


function getzero(H::Hops)
    zero_idx = zerokey(H)
    H.hoppingmatrices[:, :, zero_idx]
end

function setzero!(H::Hops, M::AbstractMatrix)
    zero_idx = zerokey(H)
    H.hoppingmatrices[:, :, zero_idx] = M
    H
end

# Hop dimension
hopdim(hops::Hops) = size(hops.hoppingmatrices, 1)

# Function to check if two bases are the same
function same_basis(h1::Hops, h2::Hops)
    return h1.basis.data == h2.basis.data
end

# Addition and subtraction with basis check
function Base.:+(h1::Hops, h2::Hops)
    if !same_basis(h1, h2)
        error("Bases of Hops instances must be the same for addition")
    end
    new_matrices = h1.hoppingmatrices .+ h2.hoppingmatrices
    Hops(h1.basis, new_matrices)
end

function Base.:-(h1::Hops, h2::Hops)
    if !same_basis(h1, h2)
        error("Bases of Hops instances must be the same for subtraction")
    end
    new_matrices = h1.hoppingmatrices .- h2.hoppingmatrices
    Hops(h1.basis, new_matrices)
end

# Scalar and matrix multiplication
function Base.:*(h::Hops, s::Number)
    new_matrices = h.hoppingmatrices .* s
    Hops(h.basis, new_matrices)
end

function Base.:*(s::Number, h::Hops)
    h * s
end

# Hops-Hops multiplication (element-wise)
function multiplyhops(h1::Hops, h2::Hops)
    if !same_basis(h1, h2)
        error("Bases of Hops instances must be the same for multiplication")
    end
    new_matrices = h1.hoppingmatrices .* h2.hoppingmatrices
    Hops(h1.basis, new_matrices)
end

# Hops-Matrix and Matrix-Hops multiplication
function multiplyhops(h1::Hops, h2::AbstractMatrix)
    new_matrices = h1.hoppingmatrices .* h2
    Hops(h1.basis, new_matrices)
end

function multiplyhops(h1::AbstractMatrix, h2::Hops)
    multiplyhops(h2, h1)
end

# In-place Hops-Hops multiplication (element-wise)
function multiplyhops!(h1::Hops, h2::Hops)
    if !same_basis(h1, h2)
        error("Bases of Hops instances must be the same for in-place multiplication")
    end
    h1.hoppingmatrices .*= h2.hoppingmatrices
    return h1
end

# In-place Hops-Scalar multiplication
function multiplyhops!(h::Hops, s::Number)
    h.hoppingmatrices .*= s
    return h
end

# Improved directsum function with basis check
function directsum(h1::Hops{K,T}, h2::Hops{K,T}) where {K,T<:AbstractHoppingMatrices}
    if !same_basis(h1, h2)
        error("Bases of Hops instances must be the same for direct sum")
    end

    d1, d2 = hopdim(h1), hopdim(h2)
    D = d1 + d2

    new_matrices = similar(h1.hoppingmatrices, D, D, size(h1.hoppingmatrices, 3))

    for i in 1:size(h1.hoppingmatrices, 3)
        new_matrices[1:d1, 1:d1, i] .= h1.hoppingmatrices[:, :, i]
        new_matrices[d1+1:end, d1+1:end, i] .= h2.hoppingmatrices[:, :, i]
    end

    Hops(h1.basis, new_matrices)
end

# Improved Kronecker product functions
function Base.kron(a, b::Hops{K,T}) where {K,T<:AbstractHoppingMatrices}
    new_matrices = similar(b.hoppingmatrices, size(a, 1) * size(b.hoppingmatrices, 1), size(a, 2) * size(b.hoppingmatrices, 2), size(b.hoppingmatrices, 3))

    for i in 1:size(b.hoppingmatrices, 3)
        new_matrices[:, :, i] .= kron(a, b.hoppingmatrices[:, :, i])
    end

    Hops(b.basis, new_matrices)
end

function Base.kron(a::Hops{K,T}, b) where {K,T<:AbstractHoppingMatrices}
    new_matrices = similar(a.hoppingmatrices, size(a.hoppingmatrices, 1) * size(b, 1), size(a.hoppingmatrices, 2) * size(b, 2), size(a.hoppingmatrices, 3))

    for i in 1:size(a.hoppingmatrices, 3)
        new_matrices[:, :, i] .= kron(a.hoppingmatrices[:, :, i], b)
    end

    Hops(a.basis, new_matrices)
end

import ..Utils

function addspin(hoppings, mode=:nospin) #::AbstractHops
    if mode == :nospin || mode == :id
        return hoppings
    elseif mode == :spinhalf || mode == :σ0
        hoppings = kron(hoppings, Utils.σ0)
    elseif mode == :σx
        hoppings = kron(hoppings, Utils.σX)
    else
        error("Do not recognize mode '$mode' in addspin(...).")
    end

    hoppings
end

decidetype(hops::Hops) = decidetype(hopdim(hops))

function decidetype(N::Int)
    if N < MAX_DENSE + 1
        return :dense
    else
        return :sparse
    end
end


function ensuretype(hops::Hops, format=:auto)
    if format == :auto
        # format = decidetype(hops) # old behaviour
        return hops
    end

    if format == :dense
        hops = DenseHops(hops)
    elseif format == :sparse
        hops = SparseHops(hops)
    end

    hops
end


function trim!(ρ::Dict{K,T}; kwargs...) where {K, T<:AbstractMatrix{ComplexF64}}
    for δL in keys(ρ)
        if all(isapprox.(ρ[δL], 0; kwargs...))
            delete!(ρ, δL)
        end
    end
    ρ
end
function trim!(ρ::Dict{K,T}; kwargs...) where {K, T<:SparseMatrixCSC{ComplexF64,Int64}}
    for δL in keys(ρ)
        ρ[δL] = SparseArrays.dropzeros!(ρ[δL])
        if all(isapprox.(ρ[δL], 0; kwargs...))
            delete!(ρ.data, δL)
        end
    end
    ρ
end

# Trimming for dense matrices
function trim!(ρ::Hops{K,<:AbstractMatrix{ComplexF64}}; kwargs...) where {K}
    # for i in 1:size(ρ.hoppingmatrices, 3)
    #     ρ.hoppingmatrices[:, :, i] .= isapprox.(ρ.hoppingmatrices[:, :, i], 0; kwargs...) .? 0 : ρ.hoppingmatrices[:, :, i]
    # end
    ρ
end

# Trimming for sparse matrices
function trim!(ρ::SparseHops{K}; kwargs...) where {K}
    for i in 1:size(ρ.hoppingmatrices, 3)
        ρ.hoppingmatrices[:, :, i] .= isapprox.(ρ.hoppingmatrices[:, :, i], 0; kwargs...) .? 0 .: ρ.hoppingmatrices[:, :, i]
        ρ.hoppingmatrices[:, :, i] = SparseArrays.dropzeros!(ρ.hoppingmatrices[:, :, i])
    end
    ρ
end

function efficientformat(ρ::Dict{K,T}) where {K,T<:AbstractMatrix{ComplexF64}}
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

function efficientzero(ρ::Dict{K,T}) where {K,T<:AbstractMatrix{ComplexF64}}
    L = length(ρ)
    @assert L > 0 "Must have at least one hopping element."

    dims = size(first(values(ρ)))
    A = zeros(ComplexF64, dims..., L) #zeros(eltype(valtype(ρ)), dims..., L)

    A, collect(keys(ρ))
end

function flexibleformat(A::AbstractArray, keylist::AbstractVector)
    Dict(L=>m for (L,m)=zip(keylist, eachslice(A; dims=3)))
end

function flexibleformat!(ρ::Dict{K,T}, A::AbstractArray, keylist::AbstractVector) where {K,T<:AbstractMatrix{ComplexF64}}
    for (j_,L)=enumerate(keylist)
        ρ[L] .= A[:,:,j_]
    end
    ρ
end