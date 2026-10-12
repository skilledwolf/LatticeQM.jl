using Test
using LinearAlgebra
using LatticeQM
using Random: Xoshiro
using LatticeQM.Spectrum: statesgrid1D, Wilson1D, wilsonweights, sectorwilson,
                          WilsonWeights1D, SectorWilson1D, statesgrid1Dslice

# Two-leg ladder of SSH chains, basis (A1, B1, A2, B2), rung tperp on A1-B2,
# bias ±V. Leg ℓ has intracell amplitude tp[ℓ] and intercell amplitude t = 1.
function ladder(tp1, tp2, tperp, V)
    h = DenseHops()
    h[[0]] = ComplexF64[V tp1 0 0; tp1 V 0 0; 0 0 -V tp2; 0 0 tp2 -V]
    h[[1]] = zeros(ComplexF64, 4, 4)
    h[[1]][1, 2] = 1.0         # A1(n) -> B1(n+1)
    h[[1]][3, 4] = 1.0         # A2(n) -> B2(n+1)
    h[[1]][1, 4] = tperp       # A1(n) -> B2(n+1) (Bernal rung)
    h[[-1]] = h[[1]]'
    h
end
τz = Diagonal([1.0, 1.0, -1.0, -1.0])

@testset "wilsonweights and sectorwilson: decoupled legs" begin
    # leg 1 topological (|tp| < t), leg 2 trivial, no rung, no bias
    h = ladder(0.4, 1.6, 0.0, 0.0)
    centres, w = WilsonWeights1D(h, 201, Matrix(τz), 1:2)
    @test sort(round.(mod.(centres .+ 0.25, 1) .- 0.25; digits=6)) == [0.0, 0.5]
    # the half-integer centre lives on leg 1 (τz = +1)
    i = argmax(abs.(mod.(centres .+ 0.25, 1) .- 0.25))
    @test isapprox(w[i], 1.0; atol=1e-6)
    p, gap = SectorWilson1D(h, 201, Matrix(τz), 1:2)
    @test isapprox(gap, 1.0; atol=1e-9)
    # Positive sector carries 1/2; negative sector carries 0.
    @test cis(2π * only(p.positive)) ≈ -1 atol=1e-12
    @test only(p.negative) ≈ 0 atol=1e-12
    # with decoupled legs the two constructions agree
    @test isapprox(sum(w .* centres) - (0.5 * 1.0), 0.0; atol=1e-6) || isapprox(abs(sum(w .* centres)), 0.5; atol=1e-6)
end

@testset "sectorwilson: coupled legs and grid refinement" begin
    h = ladder(0.4, 1.6, 0.3, 0.15)
    p, gap = SectorWilson1D(h, 301, Matrix(τz), 1:2)
    @test 0 < gap <= 1
    # phases are mod 1 and change by less than 1e-3 when the grid is refined
    p2, _ = SectorWilson1D(h, 601, Matrix(τz), 1:2)
    for sector in (:negative, :positive)
        @test maximum(abs.(mod.(p[sector] .- p2[sector] .+ 0.5, 1) .- 0.5)) < 1e-3
    end
    # the holonomy-first weights lie in [-1, 1] and the centres are mod 1
    centres, w = WilsonWeights1D(h, 301, Matrix(τz), 1:2)
    @test all(abs.(w) .<= 1 + 1e-12)
    @test all(-0.5 .<= centres .<= 0.5)
end

@testset "statesgrid1Dslice: 2D Hamiltonian at fixed momentum" begin
    lat = Geometries.honeycomb()
    H = Operators.graphene(lat)
    _, sg = statesgrid1Dslice(H, 101, [1], 1, 0.2)   # loop along k1 at k2 = 0.2
    @test size(sg) == (101, 2, 1)
    _, sg2 = statesgrid1Dslice(H, 101, [1], 2, 0.2)  # loop along k2 at k1 = 0.2
    @test size(sg2) == (101, 2, 1)
    # Zak phase of the lower graphene band along k2 at fixed k1: 0 or 1/2 mod 1
    # (chiral symmetry); away from the Dirac projections it is well defined.
    z = Wilson1D(sg[:, :, 1])
    @test isapprox(mod(2z, 1), 0.0; atol=1e-6) || isapprox(mod(2z, 1), 1.0; atol=1e-6)
end

# Change only the occupied basis at each momentum, preserving its projector.
function rotate_wilson_frames(states; seed=42)
    rng = Xoshiro(seed)
    rotated = similar(states)
    nb = size(states, 3)
    for j in axes(states, 1)
        Q = Matrix(qr(randn(rng, ComplexF64, nb, nb)).Q)
        rotated[j, :, :] = states[j, :, :] * Q
    end
    rotated
end

@testset "sectorwilson: same-sign degeneracies are subspaces" begin
    # Constant rank-two sign sectors must have zero Wilson phases, even when
    # the eigensolver is free to rotate their internal bases at every point.
    N = 101
    states = repeat(reshape(Matrix{ComplexF64}(I, 4, 4), 1, 4, 4), N, 1, 1)
    O = Diagonal([-1.0, -1.0, 1.0, 1.0])
    rotated = rotate_wilson_frames(states)
    for frames in (states, rotated)
        p, gap = sectorwilson(frames, O)
        @test p.negative ≈ zeros(2) atol=1e-12
        @test p.positive ≈ zeros(2) atol=1e-12
        @test gap ≈ 1 atol=1e-12
    end

    # A nontrivial rank-two negative sector: its Wilson eigenvalues survive
    # arbitrary occupied-space mixing despite exact degeneracy of P O P.
    H1(k) = ComplexF64[0.7 0.4+exp(-2π*im*k[1]); 0.4+exp(2π*im*k[1]) -0.7]
    H2(k) = ComplexF64[1.1 0.4+exp(-2π*im*k[1]); 0.4+exp(2π*im*k[1]) -1.1]
    _, s1 = statesgrid1D(H1, N, [1]; multimode=:serial)
    _, s2 = statesgrid1D(H2, N, [1]; multimode=:serial)
    states = zeros(ComplexF64, N, 4, 2)
    states[:, 1:2, 1] = s1[:, :, 1]
    states[:, 3:4, 2] = s2[:, :, 1]
    expected = sort([Wilson1D(s1[:, :, 1]), Wilson1D(s2[:, :, 1])])
    for frames in (states, rotate_wilson_frames(states)), sign in (-1, 1)
        p, gap = sectorwilson(frames, sign * Matrix{Float64}(I, 4, 4))
        filled, empty = sign < 0 ? (p.negative, p.positive) : (p.positive, p.negative)
        @test filled ≈ expected atol=1e-12
        @test isempty(empty)
        @test gap ≈ 1 atol=1e-12
    end
end

@testset "wilsonweights: gauge invariance and phase branch cut" begin
    # The π band can appear at either end of a sorted phase list. Transport
    # its eigenvector so its operator weight cannot swap with the zero band.
    h = ladder(0.4, 1.6, 0.0, 0.0)
    _, states = statesgrid1D(h, 101, 1:2; multimode=:serial)
    for frames in (states, rotate_wilson_frames(states))
        centres, weights = wilsonweights(frames, τz)
        iπ = argmin(abs.(cis.(2π .* centres) .+ 1))
        @test cis(2π * centres[iπ]) ≈ -1 atol=1e-12
        @test weights[iπ] ≈ 1 atol=1e-12
        @test cis(2π * centres[3-iπ]) ≈ 1 atol=1e-12
        @test weights[3-iπ] ≈ -1 atol=1e-12
    end

    h = ladder(0.4, 1.6, 0.3, 0.15)
    _, states = statesgrid1D(h, 101, 1:2; multimode=:serial)
    c, w = wilsonweights(states, τz)
    rotated = rotate_wilson_frames(states; seed=5)
    c2, w2 = wilsonweights(rotated, τz)
    @test cis.(2π .* c2) ≈ cis.(2π .* c) atol=1e-12
    @test w2 ≈ w atol=1e-12
    @test c ≈ first(Wilson1D(states)) atol=1e-12
    p, gap = sectorwilson(states, τz)
    p2, gap2 = sectorwilson(rotated, τz)
    @test cis.(2π .* p2.negative) ≈ cis.(2π .* p.negative) atol=1e-12
    @test cis.(2π .* p2.positive) ≈ cis.(2π .* p.positive) atol=1e-12
    @test gap2 ≈ gap atol=1e-12
end

@testset "Wilson loops: long products do not underflow" begin
    N, winding = 4001, 1300
    states = zeros(ComplexF64, N, 2, 1)
    for j in 1:N
        states[j, :, 1] = [sqrt(3)/2, cis(2π*winding*(j-1)/(N-1))/2]
    end
    # Each non-closing overlap equals z, and the closing overlap equals 1.
    # Its raw product underflows, but its phase is analytically well defined.
    z = (3 + cis(2π*winding/(N-1))) / 4
    expected = cis((N-1) * angle(z))
    @test abs(z)^(N-1) < floatmin(Float64)
    @test cis(2π * Wilson1D(states[:, :, 1])) ≈ expected atol=1e-10
    @test cis(2π * only(first(Wilson1D(states)))) ≈ expected atol=1e-10
    centres, weights = wilsonweights(states, Matrix{Float64}(I, 2, 2))
    @test cis(2π * only(centres)) ≈ expected atol=1e-10
    @test only(weights) ≈ 1 atol=1e-10
    p, gap = sectorwilson(states, Matrix{Float64}(I, 2, 2))
    @test cis(2π * only(p.positive)) ≈ expected atol=1e-10
    @test isempty(p.negative)
    @test gap ≈ 1 atol=1e-12
end

@testset "Operator Wilson loops: undefined inputs are rejected" begin
    # A two-dimensional constant occupied space has degenerate Wilson phases.
    states = repeat(reshape(Matrix{ComplexF64}(I, 2, 2), 1, 2, 2), 3, 1, 1)
    O = Diagonal([-1.0, 1.0])
    @test_throws r"degenerate Wilson phases" wilsonweights(states, O)
    @test_throws r"degenerate Wilson phases" wilsonweights(rotate_wilson_frames(states), O)
    @test_throws r"operator gap" sectorwilson(states, zeros(2, 2))
    @test_throws r"operator gap" sectorwilson(states, Diagonal([-1.0, 1e-16]))
    for scale in (1e-20, 1e20)
        p, gap = sectorwilson(states, scale * O)
        @test only(p.negative) ≈ 0 atol=1e-12
        @test only(p.positive) ≈ 0 atol=1e-12
        @test gap ≈ scale
    end

    # Sign changes between samples must not silently exchange sector identity.
    crossing = zeros(ComplexF64, 2, 2, 1)
    crossing[1, :, 1] = [1, 0]
    crossing[2, :, 1] = [1/2, sqrt(3)/2]
    @test_throws r"sector dimensions change" sectorwilson(crossing, O)
    crossing[2, :, 1] = [1, 1] / sqrt(2)
    @test_throws r"operator gap" sectorwilson(crossing, O)

    # An open operator gap does not make orthogonal adjacent frames usable.
    crossing[2, :, 1] = [0, 1]
    @test_throws r"singular overlap" Wilson1D(crossing)
    @test_throws r"singular overlap" Wilson1D(crossing[:, :, 1])
    @test_throws r"singular overlap" sectorwilson(crossing, Matrix{Float64}(I, 2, 2))
    @test_throws r"singular overlap" wilsonweights(crossing, Matrix{Float64}(I, 2, 2))
    for f in (sectorwilson, wilsonweights)
        @test_throws ArgumentError f(states, [1.0 0.5; 0.0 -1.0])
        @test_throws DimensionMismatch f(states, ones(3, 3))
    end
end

@testset "statesgrid1Dslice: mapping and Hamiltonian validation" begin
    H(k) = ComplexF64[0.3+k[2] 0.4+exp(-2π*im*k[1]); 0.4+exp(2π*im*k[1]) -0.3-k[2]]
    for dim in (1, 2)
        grid, states = statesgrid1Dslice(H, 13, [1], dim, 0.2; multimode=:serial)
        for (j, x) in enumerate(grid)
            k = dim == 1 ? [x, 0.2] : [0.2, x]
            v = eigen(H(k)).vectors[:, 1]
            @test states[j, :, 1] * states[j, :, 1]' ≈ v * v' atol=1e-12
        end
    end
    @test_throws ArgumentError statesgrid1Dslice(H, 5, [1], 3, 0.2)
    @test_throws ArgumentError statesgrid1Dslice(k -> [1.0 0.5; 0.0 2.0], 5, [1], 1, 0.2; multimode=:serial)
end
