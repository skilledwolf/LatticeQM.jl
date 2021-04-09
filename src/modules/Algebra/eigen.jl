###################################################################################################
# Switch between dense/sparse (see below)
###################################################################################################

geteigvals(h; format=:dense, num_bands=nothing, kwargs...) = (format==:sparse) ?  eigvals_sparse(h; num_bands=num_bands, kwargs...) : eigvals_dense(h; kwargs...)
geteigvecs(h; format=:dense, num_bands=nothing, kwargs...) = (format==:sparse) ?  eigvecs_sparse(h; num_bands=num_bands, kwargs...) : eigvecs_dense(h; kwargs...)
geteigen(h; format=:dense,   num_bands=nothing, kwargs...) = (format==:sparse) ?  eigen_sparse(h; num_bands=num_bands, kwargs...) : eigen_dense(h; kwargs...)

###################################################################################################
# Sparse eigensolver
###################################################################################################

include("eigen_sparse.jl")

# ## This is the well-tested implementation with Arpack. Unfortunately not multi-thread safe!
# eigmax_sparse(H::AbstractMatrix; kwargs...) = eigs(Hermitian(H); nev=1, which=:LR, kwargs...)[1][1] |> real
# eigmin_sparse(H::AbstractMatrix; kwargs...) = eigs(Hermitian(H); nev=1, which=:SR, kwargs...)[1][1] |> real

# eigen_sparse(M::AbstractMatrix; num_bands::Int, sigma::Float64=1e-8, which=:LM, kwargs...) = eigs(M; nev=num_bands, sigma=sigma, which=which, kwargs...)
# eigvals_sparse(args...; kwargs...) = real.((eigen_sparse(args...; kwargs...))[1])
# eigvecs_sparse(args...; kwargs...) = (eigen_sparse(args...; kwargs...))[2]

###################################################################################################
# Dense eigensolver
###################################################################################################

eigen_dense(H::AbstractMatrix, args...; kwargs...) = eigen(Matrix(H), args...; kwargs...)
eigvals_dense(H::AbstractMatrix, args...; kwargs...) = eigvals(Matrix(H), args...; kwargs...)
eigvecs_dense(H::AbstractMatrix, args...; kwargs...) = eigvecs(Matrix(H), args...; kwargs...)


###################################################################################################
# Interfaces for matrix functions
###################################################################################################

# Dense methods
eigvals_dense(h::Function; kwargs...) = k -> eigvals(Matrix(h(k)); kwargs...)
eigvecs_dense(h::Function; kwargs...) = k -> eigvecs(Matrix(h(k)); kwargs...)
f(x) = (x.values, x.vectors)
eigen_dense(h::Function; kwargs...) = k -> f(eigen(Hermitian(Matrix(h(k))); kwargs...))

# Sparse methods
eigvals_sparse(h::Function; kwargs...) = k -> eigvals_sparse(h(k); kwargs...)
eigvecs_sparse(h::Function; kwargs...) = k -> eigvecs_sparse(h(k); kwargs...)
eigen_sparse(h::Function; kwargs...) = k -> eigen_sparse(h(k); kwargs...)



