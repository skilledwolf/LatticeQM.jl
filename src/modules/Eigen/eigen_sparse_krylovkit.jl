import KrylovKit
import Random
import LinearAlgebra: I, lu, Hermitian
import SparseArrays: SparseMatrixCSC

const KRYK_TOL = 1e-12
const KRYK_MAXITER = 200

# BlockLanczos restarts keep roughly 3/5 of krylovdim, rounded down to a
# whole block. Four blocks retain at least two blocks at a restart; the old
# scalar-Lanczos sizing retained only the requested vectors and could lose
# the entire residual block before the convergence check.
_kryk_dim(num_bands::Int) = max(4*num_bands, 20)

sanitize_hermitian_sparse(H::AbstractMatrix) = H
sanitize_hermitian_sparse(H::Hermitian{T1,SparseMatrixCSC{T1,T2}}) where {T1,T2} = SparseMatrixCSC(H)

"""
    eigen_sparse(M; num_bands, sigma=1e-8, which=:LM, tol=$KRYK_TOL, kwargs...)

Hermitian sparse eigensolver via KrylovKit block Lanczos. Non-Hermitian matrices
are rejected; only roundoff-level asymmetry is accepted.

- `which=:LM` (default): eigenvalues closest to `sigma` via shift-invert.
  This is the moiré workflow — `sigma` defaults to ~0 to find low-energy bands.
- `which=:SR` / `:LR`: smallest / largest real eigenvalue, no shift-invert.

Returns `(eigvals::Vector, eigvecs::Matrix)` with energies in ascending order.
"""
function eigen_sparse(M::AbstractMatrix; num_bands::Int,
                     sigma::Float64=1e-8, which::Symbol=:LM,
                     tol::Float64=KRYK_TOL, maxiter::Int=KRYK_MAXITER, kwargs...)
    M = sanitize_hermitian_sparse(checked_hermitian(M))
    n = size(M, 1)
    1 <= num_bands <= n || throw(ArgumentError("num_bands must lie in 1:$n"))
    krylovdim = min(_kryk_dim(num_bands), n)
    alg = KrylovKit.BlockLanczos(; krylovdim, tol, maxiter)

    # A single starting vector cannot span a degenerate eigenspace. A block
    # with num_bands independent vectors resolves every multiplicity needed
    # by this request. Use a local RNG for repeatable, thread-safe results.
    rng = Random.Xoshiro(0x1a771ce9)
    x0 = KrylovKit.Block([randn(rng, float(eltype(M)), n) for _ in 1:num_bands])

    if which === :LM
        # Shift-invert: target eigenvalues of M closest to `sigma`.
        F = lu(M - sigma * I)
        op = x -> F \ x
        vals, vecs, info = KrylovKit.eigsolve(op, x0, num_bands, :LM, alg)
        info.converged < num_bands && @warn "eigen_sparse(KrylovKit): only $(info.converged)/$num_bands converged (sigma=$sigma, krylovdim=$krylovdim)"
        nk = min(num_bands, length(vals))
        nk < num_bands && @warn "eigen_sparse(KrylovKit): Krylov space terminated early — returning $nk < $num_bands eigenpairs"
        # Back-transform: eigval λ̃ of (M - σI)^{-1} → λ = σ + 1/λ̃
        realvals = sigma .+ 1 ./ real.(vals[1:nk])
        U = reduce(hcat, vecs[1:nk])
    elseif which === :SR || which === :LR
        vals, vecs, info = KrylovKit.eigsolve(M, x0, num_bands, which, alg)
        info.converged < num_bands && @warn "eigen_sparse(KrylovKit): only $(info.converged)/$num_bands converged (which=$which)"
        nk = min(num_bands, length(vals))
        nk < num_bands && @warn "eigen_sparse(KrylovKit): Krylov space terminated early — returning $nk < $num_bands eigenpairs"
        realvals = real.(vals[1:nk])
        U = reduce(hcat, vecs[1:nk])
    else
        error("eigen_sparse(KrylovKit): unsupported `which=$which` (use :LM, :SR, or :LR)")
    end

    # KrylovKit orders :LM results by distance from sigma (after the
    # shift-invert back-transform they are NOT energy-ordered) and :LR
    # descending. The dense LAPACK path returns ascending energies, and band
    # bookkeeping downstream (insertbands!, bandgap_energy, statesgrid band
    # indices) assumes that order — unsorted output scrambled band identity
    # for every `format=:sparse` band structure. Sort ascending always.
    p = sortperm(realvals)
    return realvals[p], U[:, p]
end

eigmax_sparse(H::AbstractMatrix; kwargs...) = eigen_sparse(H; num_bands=1, which=:LR, kwargs...)[1][1]
eigmin_sparse(H::AbstractMatrix; kwargs...) = eigen_sparse(H; num_bands=1, which=:SR, kwargs...)[1][1]

eigvals_sparse(args...; kwargs...) = eigen_sparse(args...; kwargs...)[1]
eigvecs_sparse(args...; kwargs...) = eigen_sparse(args...; kwargs...)[2]
