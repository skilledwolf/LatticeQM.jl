import Arpack
import SparseArrays
import LinearAlgebra: Hermitian
import SparseArrays: SparseMatrixCSC

sanitize_hermitian_sparse(H::AbstractMatrix) = H
sanitize_hermitian_sparse(H::Hermitian{T1,SparseMatrixCSC{T1,T2}}) where {T1,T2} = SparseMatrixCSC(H)

eigmax_sparse(H; kwargs...) = Arpack.eigs(sanitize_hermitian_sparse(H); nev=1, which=:LR, kwargs...)[1][1]
eigmin_sparse(H; kwargs...) = Arpack.eigs(sanitize_hermitian_sparse(H); nev=1, which=:SR, kwargs...)[1][1]

eigen_sparse(M; num_bands::Int, sigma::Float64=1e-8, which=:LM, kwargs...) = Arpack.eigs(sanitize_hermitian_sparse(M); nev=num_bands, sigma=sigma, which=which, kwargs...)
eigvals_sparse(args...; kwargs...) = (eigen_sparse(args...; kwargs...))[1]
eigvecs_sparse(args...; kwargs...) = (eigen_sparse(args...; kwargs...))[2]

# It turns out that Hermitian{...,SparseMatrixCSC{...}} causes huge problems for Arpack.eigs, so we need to convert to SparseMatrixCSC{...} first
# eigen_sparse(H::Hermitian{T1,SparseArrays.SparseMatrixCSC{T1,T2}}; kwargs...) where {T1,T2} = eigen_sparse(SparseArrays.sparse(H); kwargs...)
# eigmax_sparse(H::Hermitian{T1,SparseArrays.SparseMatrixCSC{T1,T2}}; kwargs...) where {T1,T2} = eigmax_sparse(SparseArrays.sparse(H); kwargs...)
# eigmin_sparse(H::Hermitian{T1,SparseArrays.SparseMatrixCSC{T1,T2}}; kwargs...) where {T1,T2} = eigmin_sparse(SparseArrays.sparse(H); kwargs...)