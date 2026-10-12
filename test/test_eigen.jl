using Test
using LatticeQM
using LinearAlgebra
using Random: Xoshiro
import SparseArrays: SparseMatrixCSC, sparse

@testset "Eigen: dense solves preserve the supplied operator" begin
    E = LatticeQM.Eigen
    # Include a non-normal matrix with real eigenvalues and a unitary matrix
    # close enough to I that an approximate-Hermitian heuristic would erase it.
    for A in (ComplexF64[im 0; 0 -im], [1.0 0.7; 0.0 2.0],
              Matrix(Diagonal(exp.(im .* [1e-12, -2e-12]))))
        expected = eigvals(A)
        for storage in (copy(A), sparse(A), view(copy(A), :, :))
            original = copy(storage)
            for solve in (E.geteigvals, E.geteigvals!)
                @test solve(copy(storage)) ≈ expected atol=1e-15 rtol=1e-15
            end
            for solve in (E.geteigen, E.geteigen!)
                vals, U = solve(copy(storage))
                @test vals ≈ expected atol=1e-15 rtol=1e-15
                @test A * U ≈ U * Diagonal(vals) atol=1e-15 rtol=1e-15
            end
            U = E.geteigvecs(storage)
            @test A * U ≈ U * Diagonal(expected) atol=1e-15 rtol=1e-15
            E.geteigen(storage)
            E.geteigvals(storage)
            @test storage == original
        end
    end
    # An explicit wrapper specifies the authoritative triangle, even if the
    # unused entries disagree. Conversion to dense must preserve that choice.
    A = ComplexF64[1 99; 0.3+0.2im 2]
    for storage in (A, sparse(A), view(A, :, :))
        H = Hermitian(storage, :L)
        vals, U = E.geteigen(H)
        @test H * U ≈ U * Diagonal(vals)
        @test U' * U ≈ I
    end
end

@testset "Eigen: sparse Lanczos rejects non-Hermitian operators" begin
    for A in (ComplexF64[im 0; 0 -im], [1.0 0.7; 0.0 2.0]), which in (:LM, :SR, :LR)
        @test_throws ArgumentError LatticeQM.Eigen.geteigen(sparse(A);
            format=:sparse, num_bands=1, which=which)
    end
end

@testset "Hamiltonian callers: validate before diagonalizing" begin
    E = LatticeQM.Eigen
    # A tiny absolute scale must not hide a large relative asymmetry.
    for scale in (1.0, 1e-20), A in ([1.0 0.7; 0.0 2.0], ComplexF64[im 0; 0 -im])
        H = k -> scale .* A
        ks = [0.0 0.1; 0.0 0.2]
        options = (; multimode=:serial, hidebar=true)
        @test_throws ArgumentError Spectrum.getbands(H, ks; options...)
        @test_throws ArgumentError Spectrum.getbands(H, ks, [k -> Matrix{Float64}(I, 2, 2)]; options...)
        @test_throws ArgumentError Spectrum.statesgrid1D(H, 3; multimode=:serial)
        @test_throws ArgumentError Spectrum.statesgrid(H, 3, 3; multimode=:serial)
        @test_throws ArgumentError Spectrum.berryalongpath(H, ks; multimode=:serial)
        @test_throws ArgumentError Spectrum.density(H, ks; options...)
        @test_throws ArgumentError Spectrum.ldos(H, ks, [0.0]; format=:dense, options...)
        @test_throws ArgumentError Spectrum.getdos(H, [0.0], ks; Γ=0.1, options...)
        ρ = TightBinding.Hops(Dict([0, 0] => zeros(ComplexF64, 2, 2)))
        @test_throws ArgumentError Operators.getdensitymatrix!(ρ, H, ks; options...)
        J = k -> Matrix{ComplexF64}(I, 2, 2)
        @test_throws ArgumentError LinearResponse.opticalconductivity([0.0], H, J, J, ks; options...)
    end

    # Fourier assembly can leave roundoff in either triangle or on the diagonal.
    # Keep the Hermitian solver in this case, including for sparse Lanczos.
    for T in (Float32, Float64), scale in T.((1.0, 1e-20))
        A = scale .* Complex{T}[1 0.3+0.2im; 0.3-0.2im 2]
        A[1, 1] += scale * eps(T) * im
        A[2, 1] += scale * eps(T)
        original = copy(A)
        for storage in (A, sparse(A), view(A, :, :))
            H = E.checked_hermitian(storage)
            @test H isa Hermitian
            @test parent(H) === storage
            vals, U = E.geteigen(H)
            @test norm(A * U - U * Diagonal(vals)) / norm(A) < 20eps(T)
            @test U' * U ≈ I
        end
        @test A == original
    end
    A = ComplexF64[1+eps()*im 0.3+0.2im; 0.3-0.2im+eps() 2]
    vals, U = E.geteigen(sparse(A); format=:sparse, num_bands=1, which=:SR)
    @test norm(A * U - U * Diagonal(vals)) < 1e-12
    bands = Spectrum.getbands(k -> A, zeros(2, 1); multimode=:serial, hidebar=true)
    @test bands.bands[:, 1] ≈ eigvals(Hermitian(A))
end

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
