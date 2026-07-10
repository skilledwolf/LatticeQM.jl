using Test
using LatticeQM
using LinearAlgebra
using Random: Xoshiro
import SparseArrays: SparseMatrixCSC, sparse

# Sparse eigensolver regression. The package switched from Arpack to KrylovKit;
# these tests are the safety net that catches a numerical regression on either
# direction (e.g. if the legacy Arpack solver is ever restored from git
# history, or if a future KrylovKit version changes Lanczos behaviour).
@testset "Eigen.eigen_sparse: small dense reference" begin
    # Build a small Hermitian sparse matrix with known structure and compare
    # KrylovKit's bottom-N eigenvalues to dense LAPACK ground truth.
    N = 60
    rng = Xoshiro(0x123456789abcdef)
    A = SparseMatrixCSC{ComplexF64}(sparse(Symmetric(rand(rng, N, N) .- 0.5)))
    A_dense = Matrix(A)
    ϵ_dense = sort(real.(LatticeQM.Eigen.eigen_dense(A_dense).values))

    num_bands = 10
    # :LM with shift-invert at sigma=0 — finds eigenvalues closest to 0
    ϵ_lm, U_lm = LatticeQM.Eigen.eigen_sparse(A; num_bands=num_bands, sigma=0.0)
    # Sort ascending in real for comparison
    ϵ_lm_sorted = sort(real.(ϵ_lm))
    # The 10 closest-to-zero dense eigenvalues
    closest_to_zero = sort(ϵ_dense, by=abs)[1:num_bands] |> sort
    @test ϵ_lm_sorted ≈ closest_to_zero atol=1e-8

    # :SR — smallest real eigenvalues
    ϵ_sr, _ = LatticeQM.Eigen.eigen_sparse(A; num_bands=num_bands, which=:SR)
    @test sort(real.(ϵ_sr)) ≈ ϵ_dense[1:num_bands] atol=1e-8

    # :LR — largest real eigenvalues
    ϵ_lr, _ = LatticeQM.Eigen.eigen_sparse(A; num_bands=num_bands, which=:LR)
    @test sort(real.(ϵ_lr); rev=true) ≈ reverse(ϵ_dense[end-num_bands+1:end]) atol=1e-8

    # eigmin/eigmax convenience wrappers
    @test real(LatticeQM.Eigen.eigmin_sparse(A)) ≈ ϵ_dense[1] atol=1e-8
    @test real(LatticeQM.Eigen.eigmax_sparse(A)) ≈ ϵ_dense[end] atol=1e-8
end

@testset "Eigen.eigen_sparse: residuals and orthonormality" begin
    # A real moiré-style Hamiltonian — Lanczos must converge with tight residuals
    # and return orthonormal eigenvectors. The numerical contract that matters
    # for downstream observables (Berry curvature, density matrix).
    lat = Geometries.honeycomb_twisted(5)
    T = Operators.graphene(lat)
    H = SparseMatrixCSC{ComplexF64}(T(LatticeQM.Structure.points(kpath(lat; num_points=4))[:, 2]))

    num_bands = 10
    ϵ, U = LatticeQM.Eigen.eigen_sparse(H; num_bands=num_bands)

    # Per-band residual ‖H ψ - ϵ ψ‖ — must be tighter than 1e-9 for any plausible
    # Lanczos run. Arpack used to deliver ~1e-7 here; KrylovKit gives ~1e-12.
    for b in 1:num_bands
        r = norm(H * U[:, b] - ϵ[b] * U[:, b])
        @test r < 1e-9
    end

    # Orthonormality. Arpack would fail this at ~5e-5; KrylovKit hits machine
    # precision. Downstream code (density matrix, projectors) depends on this.
    @test opnorm(U' * U - I) < 1e-10
end

# Regression: shift-invert (:LM) results come back from KrylovKit ordered by
# distance to sigma, and :LR descending — while the dense path is ascending.
# Band bookkeeping downstream indexes bands by position, so eigen_sparse must
# always return ascending energies.
@testset "Eigen.eigen_sparse: eigenvalues are ascending (band order)" begin
    rng2 = Xoshiro(0x5eed)
    N = 60
    A = rand(rng2, N, N) .- 0.5
    A = sparse(Hermitian(A + A'))
    for (which, kw) in ((:LM, (; sigma=0.2)), (:SR, (;)), (:LR, (;)))
        ϵ, U = LatticeQM.Eigen.eigen_sparse(A; num_bands=6, which=which, kw...)
        @test issorted(real.(ϵ))
        # eigenpairs still match: A u = λ u
        for m in 1:length(ϵ)
            @test norm(A * U[:, m] .- ϵ[m] .* U[:, m]) < 1e-8
        end
    end
    # determinism: same call twice gives identical results
    ϵ1, _ = LatticeQM.Eigen.eigen_sparse(A; num_bands=4)
    ϵ2, _ = LatticeQM.Eigen.eigen_sparse(A; num_bands=4)
    @test ϵ1 == ϵ2
end
