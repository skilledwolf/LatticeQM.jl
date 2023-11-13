import KrylovKit 
import LinearMaps: LinearMap

function eigen_sparse(M::AbstractMatrix; num_bands::Int, which=:LM, sigma=1e-8, kwargs...)

    # For some reason KrylovKit.eigsolve just always dies without throwing an error
    # I don't understnad what the problem is

    if which == :SM
        # use shift and invert strategy 
        F = factorize(M - sigma*I)
        Minv = LinearMap{eltype(M)}((y, x) -> (y .= F \ x), size(M, 1), ismutating=true)

        eigvals, eigvecs, info = KrylovKit.eigsolve(Minv, howmany=num_bands, which=:LM, KrylovKit.Lanczos; ishermitian=True, kwargs...)
        eigvals .= 1 ./ eigvals
    else      
        eigvals, eigvecs, info = KrylovKit.eigsolve(M, howmany=num_bands, which=which, KrylovKit.Lanczos; ishermitian=True, kwargs...)
    end

    if info.converged < num_bands
        @warn "eigen_sparse: not converged"
    end

    eigvals, eigvecs
end

eigmax_sparse(H::AbstractMatrix; kwargs...) = eigen_sparse(H, num_bands=1, which=:LR, kwargs...)[1][1]
eigmin_sparse(H::AbstractMatrix; kwargs...) = eigen_sparse(H, num_bands=1, which=:SR, kwargs...)[1][1]

eigvals_sparse(args...; kwargs...) = eigen_sparse(args...; kwargs...)[1]
eigvecs_sparse(args...; kwargs...) = eigen_sparse(args...; kwargs...)[2]