###################################################################################################
# Additional low-level interfaces (acting on arrays)
###################################################################################################

eigmax_sparse(H::AbstractMatrix; kwargs...) = eigs(H; nev=1, which=:LR, kwargs...)[1][1] |> real
eigmin_sparse(H::AbstractMatrix; kwargs...) = eigs(H; nev=1, which=:SR, kwargs...)[1][1] |> real

eigen_sparse(M::AbstractMatrix; nev::Int, sigma::Float64=1e-8, which=:LM, kwargs...) = eigs(M; nev=nev, sigma=sigma, which=which, kwargs...)
eigvals_sparse(M::AbstractMatrix; nev::Int, sigma::Float64=1e-8, which=:LM, kwargs...) = real.((eigen_sparse(M; nev=nev, sigma=sigma, which=which, kwargs...))[1])
eigvecs_sparse(M::AbstractMatrix; nev::Int, sigma::Float64=1e-8, which=:LM, kwargs...) = real.((eigen_sparse(M; nev=nev, sigma=sigma, which=which, kwargs...))[2])


###################################################################################################
# Interfaces for matrix functions
###################################################################################################

# Dense methods
eigvals_dense(h::Function; kwargs...) = k -> eigvals(h(k); kwargs...)
eigvecs_dense(h::Function; kwargs...) = k -> eigvecs(h(k); kwargs...)
# f(x) = (x.values, x.vectors)
# eigen_dense(h::Function; kwargs...) = k -> f(eigen(Hermitian(h(k)); kwargs...))
function eigen_dense(h::Function; kwargs...)
    function atK(k::AbstractVector)
        F = eigen(h(k); kwargs...) #Hermitian(h(k))
        F.values, F.vectors
    end
    atK
end

# Sparse methods
eigvals_sparse(h::Function; kwargs...) = k -> eigvals_sparse(h(k); kwargs...)
eigvecs_sparse(h::Function; kwargs...) = k -> eigvecs_sparse(h(k); kwargs...)
eigen_sparse(h::Function; kwargs...) = k -> eigen_sparse(h(k); kwargs...)

function LinearAlgebra.eigvals(h::Function; format=:dense, num_bands=nothing, kwargs...)
    if format==:dense
        eigvals_dense(h; kwargs...)
    elseif format==:sparse

        if num_bands==nothing
            eigvals_sparse(h; kwargs...)
        else
            eigvals_sparse(h; nev=num_bands, kwargs...)
        end
    end
end

function LinearAlgebra.eigvecs(h::Function; format=:dense, num_bands=nothing, kwargs...)
    if format==:dense
        eigvecs_dense(h; kwargs...)
    elseif format==:sparse

        if num_bands==nothing
            eigvecs_sparse(h; kwargs...)
        else
            eigvecs_sparse(h; nev=num_bands, kwargs...)
        end
    end
end

function LinearAlgebra.eigen(h::Function; format=:dense, num_bands=nothing, kwargs...)
    if format==:dense
        eigen_dense(h; kwargs...)
    elseif format==:sparse
        if num_bands==nothing
            eigen_sparse(h; kwargs...)
        else
            eigen_sparse(h; nev=num_bands, kwargs...)
        end
    end
end


###################################################################################################
# Interfaces for discretized k-path
###################################################################################################
# using RecursiveArrayTools
# matrixcollect(it) = convert(Array, VectorOfArray(collect(it)))

LinearAlgebra.eigvals(h::Function, ks::AbstractMatrix; kwargs...) = Base.Generator(eigvals(h; kwargs...), eachcol(ks))
LinearAlgebra.eigvecs(h::Function, ks::AbstractMatrix; kwargs...) = Base.Generator(eigvecs(h; kwargs...), eachcol(ks))
LinearAlgebra.eigen(h::Function, ks::AbstractMatrix; kwargs...) = Base.Generator(eigen(h; kwargs...), eachcol(ks))

