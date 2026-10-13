using Test
using LatticeQM
using LinearAlgebra
using Random: Xoshiro

@testset "Wilson loops: quantized SSH phases" begin
    # Intercell hopping is 1. With no staggered onsite term, inversion/chiral
    # symmetry pins the occupied-band phase to π (t < 1) or 0 (t > 1).
    for (t, expected_factor) in ((0.4, -1.0), (1.4, 1.0))
        H(k) = ComplexF64[0 t+exp(-2π*im*k[1]); t+exp(2π*im*k[1]) 0]
        _, states = Spectrum.statesgrid1D(H, 101, [1]; multimode=:serial)
        scalar_phase = Spectrum.Wilson1D(states[:, :, 1])
        matrix_phase = only(first(Spectrum.Wilson1D(states)))
        # Compare phase factors so the equivalent centres ±1/2 both pass.
        @test cis(2π * scalar_phase) ≈ expected_factor atol=1e-12
        @test cis(2π * matrix_phase) ≈ expected_factor atol=1e-12
    end
end

@testset "Wilson loops: nonquantized phases and occupied-space gauge" begin
    # The same SSH model with a staggered onsite term breaks inversion/chiral
    # symmetry without closing the gap, allowing a non-quantized Zak phase.
    H(k) = ComplexF64[0.7 0.4+exp(-2π*im*k[1]); 0.4+exp(2π*im*k[1]) -0.7]
    _, states = Spectrum.statesgrid1D(H, 101, [1]; multimode=:serial)
    phase = Spectrum.Wilson1D(states[:, :, 1])
    @test -0.21 < phase < -0.18
    phases, _ = Spectrum.Wilson1D(states)
    @test only(phases) ≈ phase atol=1e-12
    @test only(first(Spectrum.Wilson1D(reverse(states; dims=1)))) ≈ -phase atol=1e-12

    # Two independent loops with opposite phases, embedded in four orbitals.
    # Arbitrary U(2) mixing at each k must preserve the Wilson spectrum.
    N = size(states, 1)
    occupied = zeros(ComplexF64, N, 4, 2)
    occupied[:, 1:2, 1] = states[:, :, 1]
    occupied[:, 3:4, 2] = conj.(states[:, :, 1])
    rng = Xoshiro(42)
    for j in 1:N
        Q = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
        occupied[j, :, :] = occupied[j, :, :] * Q
    end
    phases, U = Spectrum.Wilson1D(occupied)
    @test phases ≈ [phase, -phase] atol=1e-12
    @test U' * U ≈ I atol=1e-12
    # Check the returned vectors against the original unitary holonomy.
    F = Matrix{ComplexF64}(I, 2, 2)
    for j in 1:N
        link = svd(occupied[j, :, :]' * occupied[mod1(j+1, N), :, :])
        F *= link.U * link.Vt
    end
    @test F * U ≈ U * Diagonal(exp.(2π*im .* phases)) atol=1e-12

    # The 2D slice API must propagate these phases (the other direction is flat).
    grid = repeat(reshape(occupied, 1, N, 4, 2), 3, 1, 1, 1)
    p1, _, p2, _ = Spectrum.Wilson2D(grid)
    @test p1 ≈ repeat([phase -phase], 3, 1) atol=1e-12
    @test maximum(abs, p2) < 1e-12

    # Nested Wannier loops use the Wilson eigenvectors, not just eigenvalues.
    # Separable tensor-product states have the same known phase along x and y.
    loop = states[1:6:end, :, 1]
    nested_phase = Spectrum.Wilson1D(loop)
    N = size(loop, 1)
    grid = zeros(ComplexF64, N, N, 8, 2)
    for i in 1:N, j in 1:N
        v = kron(loop[i, :], loop[j, :])
        grid[i, j, 1:4, 1] = v
        grid[i, j, 5:8, 2] = conj(v)
        Q = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
        grid[i, j, :, :] = grid[i, j, :, :] * Q
    end
    p1, p2 = Spectrum.WilsonWannier2D(grid)
    @test p1 ≈ [nested_phase, -nested_phase] atol=1e-12
    @test p2 ≈ [nested_phase, -nested_phase] atol=1e-12
end

@testset "Berry: unsupported sparse mode fails before evaluating H" begin
    H = k -> error("must not construct a dense Hamiltonian")
    @test_throws ArgumentError Spectrum.statesgrid(H, 3, 3, [1]; format=:sparse)
    @test_throws ArgumentError Spectrum.statesgrid1D(H, 3, [1]; format=:sparse)
    @test_throws ArgumentError Spectrum.berryalongpath(H, zeros(2, 3); format=:sparse)
end
