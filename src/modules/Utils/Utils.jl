"""
General functions and definitions that don't necessarily make sense anywhere else.
"""
# Note: This module must NOT depend on any other module of the package.
module Utils

    using LinearAlgebra

    # include("DummySave.jl")
    # export DummySave

    # include("legacymacros.jl")
    # export @legacyalias, @legacymoved, @legacyremoved

    include("Context.jl")
    import .Context

    include("pycall.jl")

    include("scalar2vector.jl")

    include("fermidirac.jl")
    include("paulimatrices.jl")
    export σ0, σ1, σ2, σ3, σX, σY, σZ, σUP, σDOWN, σPLUS, σMINUS, σs, spinorrotation

    # Just pass through different dense types
    import SharedArrays, SparseArrays
    dense(A::Array) = A
    dense(A::SharedArrays.SharedArray) = A
    dense(A::Hermitian{T,Matrix{T}}) where {T} = A
    dense(A::Hermitian{T,SharedArrays.SharedMatrix{T}}) where {T} = A
    dense(A::SubArray) = Array(A) # for now dense creates copies of views instead of returning the view
    # Create dense copies for sparse cases
    dense(A::SparseArrays.SparseMatrixCSC) = Array(A)
    dense(A::Hermitian{T,SparseArrays.SparseMatrixCSC{T,K}}) where {T,K} = Hermitian(Array(A))

<<<<<<< HEAD
    densecopy(A::Array) = Array(A)
    densecopy(A::SubArray) = Array(A)
=======
>>>>>>> origin/master

    export dense


    getelectronsector(H::Function) = H
    getelectronsector(H::AbstractMatrix) = H
    copyelectronsector(H::Function) = H
    copyelectronsector(H::AbstractMatrix) = copy(H)

    """
        padvec(v::AbstractVector, d::Int)

    Make sure Vector v has length d, pad with zeros if needed.
    """
    function padvec(v::AbstractVector, d::Int)
        L = length(v)
        @assert L <= d "Vector exceeds specified length $d."
        (L==d) ? v : vcat(v,zeros(d-L))
    end
end
