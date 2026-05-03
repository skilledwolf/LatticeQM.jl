module Eigen 

import LinearAlgebra
import LinearAlgebra: Hermitian
import SharedArrays
import SparseArrays

import LatticeQM.Utils: dense

###################################################################################################
# Switch between dense/sparse (see below)
###################################################################################################

# function geteigvals(h; format=:dense, kwargs...)
#     if format == :sparse
#         return let h0 = sanatize_hermitian_sparse(h)
#             eigvals_sparse(h0; kwargs...)
#         end
#     else
#         return let h0 = dense(h)
#             LinearAlgebra.eigvals(h0; kwargs...)
#         end
#     end
# end

# function geteigvecs(h; format=:dense, kwargs...)
#     if format == :sparse
#         return let h0 = sanatize_hermitian_sparse(h)
#             eigvecs_sparse(h0; kwargs...)
#         end
#     else
#         return let h0 = dense(h0)
#             LinearAlgebra.eigvecs(h0; kwargs...)
#         end
#     end
# end

# function geteigen(h; format=:dense, kwargs...)
#     if format == :sparse
#         return let h0 = sanatize_hermitian_sparse(h)
#             eigen_sparse(h0; kwargs...)
#         end
#     else
#         return let h0 = dense(h)
#             LinearAlgebra.eigen(h0; kwargs...)
#         end
#     end
# end


eigvals_dense(M; kwargs...) = LinearAlgebra.eigvals(_hermitianize(dense(M)); kwargs...)
eigvecs_dense(M; kwargs...) = LinearAlgebra.eigvecs(_hermitianize(dense(M)); kwargs...)
eigen_dense(M; kwargs...) = LinearAlgebra.eigen(_hermitianize(dense(M)); kwargs...)

# In-place dense variants. `eigvals!`/`eigen!` overwrite their input but skip
# the dense() copy that the non-bang versions allocate on every call. They're
# the right path inside hot k-loops where the caller has already filled a
# scratch Hcache via `Spectrum._build_H!`.
eigvals_dense!(M; kwargs...) = LinearAlgebra.eigvals!(_hermitianize(dense(M)); kwargs...)
eigen_dense!(M; kwargs...) = LinearAlgebra.eigen!(_hermitianize(dense(M)); kwargs...)

# Hamiltonians in this package are physically Hermitian (tight-binding, BdG,
# Floquet). Wrap dense matrices in `Hermitian(...)` so LAPACK uses the
# Hermitian solver (`heevr`) instead of the generic `geevx` path. The
# Hermitian path is faster *and* numerically stable on borderline matrices —
# `LAPACKException(3)` failures were observed on BdG sweeps with the generic
# solver. `Hermitian` reads from the upper triangle and ignores the lower,
# which is fine for our user-controlled Hamiltonians.
@inline _hermitianize(M::AbstractMatrix) = LinearAlgebra.Hermitian(M)
@inline _hermitianize(M::LinearAlgebra.Hermitian) = M

geteigvals(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvals_sparse(args...; kwargs...) : eigvals_dense(args...; kwargs...)
geteigvecs(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvecs_sparse(args...; kwargs...) : eigvecs_dense(args...; kwargs...)
geteigen(args...; format=:dense, kwargs...) = (format == :sparse) ? eigen_sparse(args...; kwargs...) : eigen_dense(args...; kwargs...)
geteigen!(args...; format=:dense, kwargs...) = (format == :sparse) ? eigen_sparse(args...; kwargs...) : eigen_dense!(args...; kwargs...)
geteigvals!(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvals_sparse(args...; kwargs...) : eigvals_dense!(args...; kwargs...)

# function geteigvals(h; format=:dense, kwargs...)
#     if format == :sparse
#         return let h0 = sanatize_hermitian_sparse(h)
#             eigvals_sparse(h0; kwargs...)
#         end
#     else
#         return let h0 = dense(h)
#             LinearAlgebra.eigvals(dense(h0); kwargs...)
#         end
#     end
# end

# function geteigvecs(h; format=:dense, kwargs...)
#     if format == :sparse
#         return let h0 = sanatize_hermitian_sparse(h)
#             eigvecs_sparse(h0; kwargs...)
#         end
#     else
#         return let h0 = dense(h0)
#             LinearAlgebra.eigvecs(dense(h0); kwargs...)
#         end
#     end
# end

# function geteigen(h::AbstractMatrix; format=:dense, kwargs...)
#     if format == :sparse
#         eigen_sparse(sanatize_hermitian_sparse(h); kwargs...)
#     else
#         LinearAlgebra.eigen(h0; kwargs...)
#     end
# end

###################################################################################################
# Sparse eigensolver
###################################################################################################

# Sparse Hermitian eigensolver. KrylovKit (default) is pure-Julia, thread-safe,
# and gives tighter residuals than Arpack — see test_eigen.jl for the regression.
# To temporarily revert to Arpack (e.g. to debug a numerical regression), swap
# the include to "eigen_sparse_arpack.jl"; both files implement the same API.
include("eigen_sparse_krylovkit.jl")
# include("eigen_sparse_arpack.jl")  # legacy fallback, NOT thread-safe
# include("eigen_sparse_julia.jl")   # ArnoldiMethod fallback (incomplete)

end # module Eigen