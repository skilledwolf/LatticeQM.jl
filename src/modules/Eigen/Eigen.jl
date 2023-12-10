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


eigvals_dense(M; kwargs...) = LinearAlgebra.eigvals(dense(M); kwargs...)
eigvecs_dense(M; kwargs...) = LinearAlgebra.eigvecs(dense(M); kwargs...)
eigen_dense(M; kwargs...) = LinearAlgebra.eigen(dense(M); kwargs...)

geteigvals(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvals_sparse(args...; kwargs...) : eigvals_dense(args...; kwargs...)
geteigvecs(args...; format=:dense, kwargs...) = (format == :sparse) ? eigvecs_sparse(args...; kwargs...) : eigvecs_dense(args...; kwargs...)
geteigen(args...; format=:dense, kwargs...) = (format == :sparse) ? eigen_sparse(args...; kwargs...) : eigen_dense(args...; kwargs...)

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

# include("eigen_sparse_krylovkit.jl")
# include("eigen_sparse_julia.jl")
include("eigen_sparse_arpack.jl") # Arpack implementation. Not multi-threading safe with julia.

end # module Eigen