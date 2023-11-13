import Arpack
eigmax_sparse(H::AbstractMatrix; kwargs...) = Arpack.eigs(H; nev=1, which=:LR, kwargs...)[1][1] #|> real
eigmin_sparse(H::AbstractMatrix; kwargs...) = Arpack.eigs(H; nev=1, which=:SR, kwargs...)[1][1] #|> real

eigen_sparse(M::AbstractMatrix; num_bands::Int, sigma::Float64=1e-8, which=:LM, kwargs...) = Arpack.eigs(M; nev=num_bands, sigma=sigma, which=which, kwargs...)
eigvals_sparse(args...; kwargs...) = (eigen_sparse(args...; kwargs...))[1]
eigvecs_sparse(args...; kwargs...) = (eigen_sparse(args...; kwargs...))[2]