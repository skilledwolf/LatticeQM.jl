module Eigen 

import LinearAlgebra
import LinearAlgebra: Hermitian
import SharedArrays
import SparseArrays

import LatticeQM.Utils: dense

###################################################################################################
# Switch between dense/sparse (see below)
###################################################################################################

# Preserve the operator: Wilson holonomies, for example, are unitary rather
# than Hermitian. Julia dispatches on explicit matrix structure and detects
# exact Hermiticity itself for ordinary dense matrices.
eigvals_dense(M; kwargs...) = LinearAlgebra.eigvals(dense(M); kwargs...)
eigvecs_dense(M; kwargs...) = LinearAlgebra.eigvecs(dense(M); kwargs...)
eigen_dense(M; kwargs...) = LinearAlgebra.eigen(dense(M); kwargs...)

# In-place dense variants. `eigvals!`/`eigen!` overwrite their input but skip
# the solver's input copy for dense matrices. They're the right path inside
# hot k-loops where the caller has already filled a
# scratch Hcache via `Spectrum._build_H!`.
eigvals_dense!(M; kwargs...) = LinearAlgebra.eigvals!(dense(M); kwargs...)
eigen_dense!(M; kwargs...) = LinearAlgebra.eigen!(dense(M); kwargs...)

"""
    checked_hermitian(M)

Internal boundary check for Hamiltonian observables and Hermitian Lanczos.
Accept roundoff from Bloch/BdG assembly (relative Frobenius-norm discrepancy
at most `100eps` of the real floating element type, with no absolute floor),
then use the upper triangle. Materially non-Hermitian inputs throw an
`ArgumentError`. An explicit `Hermitian` wrapper already defines the operator
and its authoritative triangle, so it is preserved.

This retains the Hermitian LAPACK path for approximately Hermitian BdG
matrices, where generic `geevx` previously failed with `LAPACKException(3)`.
Generic dense eigensolvers must not use this check: even tiny non-Hermitian
parts can carry meaningful Wilson phases.
"""
@inline checked_hermitian(M::Hermitian) = M
function checked_hermitian(M::AbstractMatrix)
    if size(M, 1) == size(M, 2) &&
       (LinearAlgebra.ishermitian(M) ||
        isapprox(M, M'; rtol=100eps(real(float(one(eltype(M))))), atol=0))
        return Hermitian(M)
    end
    throw(ArgumentError("expected a Hermitian Hamiltonian/operator (up to roundoff); use Eigen.geteigen(M; format=:dense) for a general operator"))
end

geteigvals(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvals_sparse(args...; kwargs...) : eigvals_dense(args...; kwargs...)
geteigvecs(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvecs_sparse(args...; kwargs...) : eigvecs_dense(args...; kwargs...)
geteigen(args...; format=:dense, kwargs...) = (format == :sparse) ? eigen_sparse(args...; kwargs...) : eigen_dense(args...; kwargs...)
geteigen!(args...; format=:dense, kwargs...) = (format == :sparse) ? eigen_sparse(args...; kwargs...) : eigen_dense!(args...; kwargs...)
geteigvals!(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvals_sparse(args...; kwargs...) : eigvals_dense!(args...; kwargs...)

###################################################################################################
# Sparse eigensolver
###################################################################################################

# Sparse Hermitian eigensolver. KrylovKit is pure-Julia, thread-safe, and gives
# tighter residuals than the Arpack solver it replaced — see test_eigen.jl for
# the regression. The legacy implementations (eigen_sparse_arpack.jl, Arpack,
# NOT thread-safe; eigen_sparse_julia.jl, ArnoldiMethod, incomplete) were
# removed together with their dependencies; recover them from git history if a
# fallback is ever needed for debugging. Both implemented the same API.
include("eigen_sparse_krylovkit.jl")

end # module Eigen
