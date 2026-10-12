"""
General functions and definitions that don't necessarily make sense anywhere else.
"""
# Note: This module must NOT depend on any other module of the package.
module Utils

    using LinearAlgebra

    include("scalar2vector.jl")

    include("fermidirac.jl")
    include("paulimatrices.jl")
    export σ0, σ1, σ2, σ3, σX, σY, σZ, σUP, σDOWN, σPLUS, σMINUS, σs, spinorrotation

    # Just pass through different dense types
    import SharedArrays
    dense(A::Array) = A
    dense(A::SharedArrays.SharedArray) = A
    dense(A::SubArray) = Array(A) # for now dense creates copies of views instead of returning the view
    # Materialize other storage, preserving an explicit Hermitian contract and
    # the selected triangle even for sparse matrices and views.
    dense(A::AbstractMatrix) = Matrix(A)
    dense(A::Hermitian) = Hermitian(dense(parent(A)), A.uplo == 'U' ? :U : :L)

    densecopy(A::Array) = Array(A)
    densecopy(A::SubArray) = Array(A)
    densecopy(A::SharedArrays.SharedArray) = Array(A)

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
