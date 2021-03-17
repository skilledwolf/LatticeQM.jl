using ArnoldiMethod, LinearAlgebra, LinearMaps

################################################################################
# Direct interfaces to eigenvalues or eigenvectors
################################################################################
eigvals_sparse(args...; kwargs...) = real.((eigen_sparse(args...; kwargs...))[1])
eigvecs_sparse(args...; kwargs...) = (eigen_sparse(args...; kwargs...))[2]


################################################################################
# Largest/smallest eigenvalue
################################################################################

function eigmax_sparse(H::AbstractMatrix; kwargs...)
    decomp, _ = partialschur(A, nev=1, tol=sigma, which=LR(), kwargs...)
    λs, _ = partialeigen(decomp)

    real(λs[1])
end

function eigmin_sparse(H::AbstractMatrix; sigma=1e-6, kwargs...)
    decomp, _ = partialschur(A, nev=1, tol=sigma, which=SR(), kwargs...)
    λs, _ = partialeigen(decomp)

    real(λs[1])
end

################################################################################
# Eigen problem with shift and invert
################################################################################
# import LinearAlgebra
# import SuiteSparse
# LinearAlgebra.ldiv!(Y, F::Union{SuiteSparse.SPQR.QRSparse,SuiteSparse.CHOLMOD.Factor}, X) = copyto!(Y, F\X) # Warning! This potentially allocates significant memory!

# Factorizes A and builds a linear map that applies inv(A) to a vector.
function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) ->(y.=F\x) , size(A,1), ismutating=true) # (y, x) ->  ldiv!(y, F, x)
end

function eigen_sparse(A::AbstractMatrix; num_bands, kwargs...)

    # Target the largest eigenvalues of the inverted problem
    decomp, _ = partialschur(construct_linear_map(A), nev=num_bands, which=LM(), kwargs...) #nev=4, tol=1e-5, restarts=100, which=LM()
    λs_inv, X = partialeigen(decomp)

    # Eigenvalues have to be inverted to find the smallest eigenvalues of the non-inverted problem.
    λs = 1 ./ λs_inv

    λs, X
end

################################################################################
# Generalized eigen problem
################################################################################

struct ShiftAndInvert{TA,TB,TT}
    A_lu::TA
    B::TB
    temp::TT
end

function (M::ShiftAndInvert)(y,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A_lu, M.temp)
end

function construct_linear_map(A,B)
    a = ShiftAndInvert(factorize(A),B,Vector{eltype(A)}(undef, size(A,1)))
    LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
end
 
function eigen_sparse(A::AbstractMatrix, B::AbstractMatrix; num_bands, kwargs...)

    # Target the largest eigenvalues of the inverted problem
    decomp, _ = partialschur(construct_linear_map(A, B), nev=num_bands, which=LM(), kwargs...)
    λs_inv, X = partialeigen(decomp)

    # Eigenvalues have to be inverted to find the smallest eigenvalues of the non-inverted problem.
    λs = 1 ./ λs_inv

    λs, X
end