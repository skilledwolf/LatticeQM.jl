# Threaded sparse-eigensolver experiment + accuracy check.
#
# Run as:  julia --project=... extra/benchmarks/bench_kryk.jl <case> <mode>
#
# <case> ∈ {tbg_n7, tbg_n11}
# <mode> ∈ {accuracy, serial_kryk, threaded_kryk}
#
# - accuracy: build H(k) at a few k-points, solve with both Arpack and
#   KrylovKit, report max abs eigenvalue difference and max abs(1 - |⟨ψ_kryk,
#   ψ_arpack⟩|) for the lowest |E| eigenpairs.
# - serial_kryk / threaded_kryk: full bandmatrix-equivalent run using KrylovKit.

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

using LatticeQM
import KrylovKit
import Arpack
import SparseArrays: SparseMatrixCSC, sparse
import LatticeQM.Spectrum: handleprojector, insertbands_bandexpvals_k!
import LatticeQM.Eigen

# Thread-safe sparse eigensolver via KrylovKit + shift-invert (target eigenvalues
# closest to `sigma`). Mirrors Eigen.geteigen output: (eigvals, eigvecs::Matrix).
function eigen_sparse_kryk(M::SparseMatrixCSC{ComplexF64}; num_bands::Int, sigma::Float64=1e-8)
    F = lu(M - sigma*I)
    x0 = randn(ComplexF64, size(M, 1))
    op = x -> F \ x
    vals, vecs, info = KrylovKit.eigsolve(op, x0, num_bands, :LM,
                                          KrylovKit.Lanczos(; krylovdim=max(2*num_bands+10, 20),
                                                              tol=1e-12))
    info.converged < num_bands && @warn "kryk: only $(info.converged) of $num_bands converged"
    realvals = sigma .+ 1 ./ vals[1:num_bands]
    U = reduce(hcat, vecs[1:num_bands])
    realvals, U
end

case = ARGS[1]
mode = ARGS[2]

n = parse(Int, replace(case, "tbg_n" => ""))
lat = Geometries.honeycomb_twisted(n)
T = Operators.graphene(lat)
norbital = LatticeQM.Structure.Lattices.countorbitals(lat)
nks = norbital > 700 ? 16 : 32
ks_path = kpath(lat; num_points=nks)
ks = LatticeQM.Structure.points(ks_path)

num_bands = 10

# Build a helper that materializes H(k) as SparseMatrixCSC{ComplexF64}.
makeH(k) = SparseMatrixCSC{ComplexF64}(T(k))

function accuracy_check(ks, num_bands)
    sample = [1, fld(size(ks, 2), 3), 2*fld(size(ks, 2), 3), size(ks, 2)]
    max_eval_diff = 0.0
    max_subspace_err = 0.0   # ||P_arpack - P_kryk||_F over the low-energy subspace
    max_residual_a = 0.0     # ||H ψ_a - ϵ_a ψ_a|| (Arpack residual)
    max_residual_k = 0.0     # ||H ψ_k - ϵ_k ψ_k|| (KrylovKit residual)
    for j in sample
        H = makeH(@view ks[:, j])
        ϵ_a, U_a = Eigen.geteigen(H; format=:sparse, num_bands=num_bands)
        ϵ_k, U_k = eigen_sparse_kryk(H; num_bands=num_bands)
        pa = sortperm(real.(ϵ_a)); pk = sortperm(real.(ϵ_k))
        ϵ_a, U_a = ϵ_a[pa], U_a[:, pa]
        ϵ_k, U_k = ϵ_k[pk], U_k[:, pk]
        max_eval_diff = max(max_eval_diff, maximum(abs.(real.(ϵ_a) .- real.(ϵ_k))))
        # Subspace projectors are basis-independent — robust to eigenvector
        # gauge and to degeneracies (high-symmetry k-points have band crossings).
        P_a = U_a * U_a'
        P_k = U_k * U_k'
        max_subspace_err = max(max_subspace_err, opnorm(P_a - P_k))
        # Residuals: how well does each solver actually solve H ψ = ϵ ψ?
        for b in 1:num_bands
            r_a = H * U_a[:, b] - ϵ_a[b] * U_a[:, b]
            r_k = H * U_k[:, b] - ϵ_k[b] * U_k[:, b]
            max_residual_a = max(max_residual_a, norm(r_a))
            max_residual_k = max(max_residual_k, norm(r_k))
        end
    end
    max_eval_diff, max_subspace_err, max_residual_a, max_residual_k
end

if mode == "accuracy"
    d, s, ra, rk = accuracy_check(ks, num_bands)
    println("$case,accuracy",
            ",eval_diff=", round(d; sigdigits=4),
            ",subspace_diff=", round(s; sigdigits=4),
            ",residual_arpack=", round(ra; sigdigits=4),
            ",residual_kryk=", round(rk; sigdigits=4))
    exit(0)
end

projectors = handleprojector()
bands = zeros(Float64, num_bands, size(ks, 2))
obs = zeros(Float64, num_bands, size(ks, 2), 0)

# warm-up
let
    H0 = makeH(@view ks[:, 1])
    eigen_sparse_kryk(H0; num_bands=num_bands)
end

GC.gc(); GC.gc()
t = @timed begin
    if mode == "serial_kryk"
        for j in 1:size(ks, 2)
            H0 = makeH(@view ks[:, j])
            ϵs, U = eigen_sparse_kryk(H0; num_bands=num_bands)
            @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], ϵs, U, ks[:, j], projectors)
        end
    elseif mode == "threaded_kryk"
        Threads.@threads for j in 1:size(ks, 2)
            H0 = makeH(@view ks[:, j])
            ϵs, U = eigen_sparse_kryk(H0; num_bands=num_bands)
            @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], ϵs, U, ks[:, j], projectors)
        end
    else
        error("unknown mode $mode")
    end
end

elapsed = t.time
allocs_mb = t.bytes / 1024^2
rss_mb = Sys.maxrss() / 1024^2

println(join((case, mode, Threads.nthreads(), 1, LinearAlgebra.BLAS.get_num_threads(),
              size(ks, 2), norbital, round(elapsed, digits=3), round(rss_mb, digits=1),
              round(allocs_mb, digits=1)), ","))
