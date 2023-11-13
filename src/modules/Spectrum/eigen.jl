###################################################################################################
# Switch between dense/sparse (see below)
###################################################################################################

geteigvals(h; format=:dense, kwargs...) = (format==:sparse) ?  eigvals_sparse(h; kwargs...) : eigvals_dense(h; kwargs...)
geteigvecs(h; format=:dense, kwargs...) = (format==:sparse) ?  eigvecs_sparse(h; kwargs...) : eigvecs_dense(h; kwargs...)
geteigen(h; format=:dense,   kwargs...) = (format==:sparse) ?  eigen_sparse(h; kwargs...) : eigen_dense(h; kwargs...)

###################################################################################################
# Sparse eigensolver
###################################################################################################

# include("eigen_sparse_krylovkit.jl")
# include("eigen_sparse_julia.jl")
include("eigen_sparse_arpack.jl") # Arpack implementation. Not multi-threading safe with julia.

###################################################################################################
# Dense eigensolver
###################################################################################################

import LinearAlgebra
import SharedArrays: sdata

eigen_dense(H::AbstractMatrix, args...; kwargs...) = (F=LinearAlgebra.eigen(H, args...; kwargs...); (F.values, F.vectors))
eigvals_dense(H::AbstractMatrix, args...; kwargs...) = LinearAlgebra.eigvals(H, args...; kwargs...)
eigvecs_dense(H::AbstractMatrix, args...; kwargs...) = LinearAlgebra.eigvecs(H, args...; kwargs...)


###################################################################################################
# Interfaces for matrix functions
###################################################################################################

# # Dense methods
# eigvals_dense(h::Function; kwargs...) = k -> eigvals(Matrix(h(k)); kwargs...)
# eigvecs_dense(h::Function; kwargs...) = k -> eigvecs(Matrix(h(k)); kwargs...)
# f(x) = (x.values, x.vectors)
# eigen_dense(h::Function; kwargs...) = k -> f(eigen(Hermitian(Matrix(h(k))); kwargs...))

# # Sparse methods
# eigvals_sparse(h::Function; kwargs...) = k -> eigvals_sparse(h(k); kwargs...)
# eigvecs_sparse(h::Function; kwargs...) = k -> eigvecs_sparse(h(k); kwargs...)
# eigen_sparse(h::Function; kwargs...) = k -> eigen_sparse(h(k); kwargs...)



