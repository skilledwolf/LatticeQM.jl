using Test
using LinearAlgebra
using LatticeQM
using Random: Xoshiro
using LatticeQM.Spectrum: statesgrid1D, Wilson1D, wilsonweights, sectorwilson,
                          WilsonWeights1D, SectorWilson1D, statesgrid1Dslice,
                          Wilson2D, WilsonWannier2D

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
    # Zak phase of the lower graphene band along k1 at fixed k2: 0 or 1/2 mod 1
    # (chiral symmetry); away from the Dirac projections it is well defined.
    z = Wilson1D(sg[:, :, 1])
    @test isapprox(mod(2z, 1), 0.0; atol=1e-6) || isapprox(mod(2z, 1), 1.0; atol=1e-6)
end

# Change only the occupied basis at each momentum, preserving its projector.
function rotate_wilson_frames(states; seed=42)
    rng = Xoshiro(seed)
    rotated = similar(states)
    nb = size(states, ndims(states))
    flat = reshape(states, :, size(states, ndims(states)-1), nb)
    out = reshape(rotated, size(flat))
    for j in axes(flat, 1)
        Q = Matrix(qr(randn(rng, ComplexF64, nb, nb)).Q)
        out[j, :, :] = flat[j, :, :] * Q
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
    states = repeat(reshape(Matrix{ComplexF64}(I, 2, 2), 1, 2, 2), 3, 1, 1)
    O = Diagonal([-1.0, 1.0])
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

@testset "Wilson weights: degenerate subspaces" begin
    # Both legs in the same SSH phase have coincident Wilson eigenvalues.
    # The layer-weight trace is zero, while the charge weight counts states.
    for (t1, t2, factor) in ((0.4, 0.6, -1), (1.4, 1.6, 1))
        h = ladder(t1, t2, 0.0, 0.0)
        kgrid, states = statesgrid1D(h, 61, 1:2; multimode=:serial)
        for frames in (states, rotate_wilson_frames(states))
            result = wilsonweights(frames, τz; kgrid)
            @test result.multiplicities == [2]
            @test cis(2π * only(result.centres)) ≈ factor atol=1e-12
            @test only(result.weights) ≈ 0 atol=1e-12
            charge = wilsonweights(frames, Matrix{Float64}(I, 4, 4); kgrid)
            @test charge.weights ≈ [2] atol=1e-12
            _, U = Wilson1D(frames)
            @test U' * U ≈ I atol=1e-12
        end
    end

    # A degenerate π pair and an isolated zero phase: preserve block identity
    # under general U(3) gauge changes, including numerical splitting at ±π.
    _, pair = statesgrid1D(ladder(0.4, 0.6, 0.0, 0.0), 61, 1:2; multimode=:serial)
    H(k) = ComplexF64[0 1.6+cis(-2π*k[1]); 1.6+cis(2π*k[1]) 0]
    _, trivial = statesgrid1D(H, 61, [1]; multimode=:serial)
    states = zeros(ComplexF64, 61, 6, 3)
    states[:, 1:4, 1:2] = pair
    states[:, 5:6, 3] = trivial[:, :, 1]
    O = Diagonal([1.0, 1.0, 2.0, 2.0, 4.0, 4.0])
    for frames in (states, rotate_wilson_frames(states; seed=31))
        result = wilsonweights(frames, O)
        iπ = argmin(abs.(cis.(2π .* result.centres) .+ 1))
        @test result.multiplicities[iπ] == 2
        @test result.weights[iπ] ≈ 3 atol=1e-12
        @test result.multiplicities[3-iπ] == 1
        @test result.weights[3-iπ] ≈ 4 atol=1e-12
    end
end

@testset "Wilson weights: quadrature and explicit occupancy" begin
    # The occupied spinor has <σz> = 0.2 + 0.4cos(2πk), whose mean is 0.2.
    function average_model(shift)
        k -> begin
            x = 2π*(k[1]+shift)
            m = 0.2 + 0.4cos(x)
            v = ComplexF64[sqrt((1+m)/2), cis(x)*sqrt((1-m)/2)]
            Matrix{ComplexF64}(I, 2, 2) - 2v*v'
        end
    end
    O = Diagonal([1.0, -1.0])
    for N in (21, 101), shift in (0.0, 0.5)
        H = average_model(shift)
        result = WilsonWeights1D(H, N, O, [1]; multimode=:serial)
        @test only(result.weights) ≈ 0.2 atol=1e-12
        kgrid, states = statesgrid1D(H, N, [1]; multimode=:serial)
        @test wilsonweights(rotate_wilson_frames(states), O; kgrid).weights ≈ result.weights atol=1e-12
        @test wilsonweights(states[1:end-1, :, :], O).weights ≈ result.weights atol=1e-12
    end
    H = average_model(0.0)
    for f in (WilsonWeights1D, SectorWilson1D)
        @test_throws MethodError f(H, 21, O)
        @test_throws r"bandindices explicitly" f(H, 21, O, Int[])
        @test_throws ArgumentError f(H, 1, O, [1])
    end
    _, states = statesgrid1D(H, 5, [1]; multimode=:serial)
    @test_throws DimensionMismatch wilsonweights(states, O; kgrid=[0, 1])
    @test_throws ArgumentError wilsonweights(states, O; kgrid=[0, 0.3, 0.2, 0.8, 1])
    @test_throws ArgumentError wilsonweights(states, O; kgrid=range(0, 2; length=5))

    # A nonuniform grid represents a piecewise-linear expectation whose
    # exact integral is 0.19. Including the repeated endpoint changes nothing.
    kgrid = [0.0, 0.1, 0.4, 0.9, 1.0]
    m = [0.6, 0.3, -0.2, 0.5, 0.6]
    for j in eachindex(kgrid)
        states[j,:,1] = [sqrt((1+m[j])/2), cis(2π*kgrid[j])*sqrt((1-m[j])/2)]
    end
    @test only(wilsonweights(states, O; kgrid).weights) ≈ 0.19 atol=1e-12
    @test only(wilsonweights(states[1:4,:,:], O; kgrid=kgrid[1:4]).weights) ≈ 0.19 atol=1e-12
end

@testset "Wilson weights and sectors: transport identities" begin
    _, states = statesgrid1D(ladder(0.4, 1.6, 0.3, 0.15), 81, 1:2; multimode=:serial)
    c, w = wilsonweights(states, τz)
    @test sum(w) ≈ sum(real(tr(states[j,:,:]' * τz * states[j,:,:])) for j=1:81)/81 atol=1e-12
    cr, wr = wilsonweights(reverse(states; dims=1), τz)
    @test cr ≈ reverse(-c) atol=1e-12
    @test wr ≈ reverse(w) atol=1e-12
    cs, ws = wilsonweights(circshift(states, (23,0,0)), τz)
    @test cs ≈ c atol=1e-12
    @test ws ≈ w atol=1e-12
    p, gap = sectorwilson(states, τz)
    pn, gapn = sectorwilson(states, -τz)
    @test p.negative ≈ pn.positive atol=1e-12
    @test p.positive ≈ pn.negative atol=1e-12
    @test gap ≈ gapn atol=1e-12
end

@testset "Nested Wilson loops: phase cuts and degenerate sectors" begin
    N = 13
    constant = zeros(ComplexF64, N, N, 4, 2)
    constant[:, :, 1, 1] .= 1
    constant[:, :, 2, 2] .= 1
    for frames in (constant, rotate_wilson_frames(constant))
        p1, p2 = WilsonWannier2D(frames)
        @test p1 ≈ [0] atol=1e-12
        @test p2 ≈ [0] atol=1e-12
    end

    ssh(t, m=0.0) = k -> ComplexF64[m t+cis(-2π*k[1]); t+cis(2π*k[1]) -m]
    _, topo = statesgrid1D(ssh(0.4), N+1, [1]; multimode=:serial)
    _, trivial = statesgrid1D(ssh(1.6), N+1, [1]; multimode=:serial)
    grid = zeros(ComplexF64, N, N, 8, 2)
    for i in 1:N, j in 1:N
        grid[i,j,1:4,1] = kron(topo[i,:,1], topo[j,:,1])
        grid[i,j,5:8,2] = kron(trivial[i,:,1], trivial[j,:,1])
    end
    for frames in (grid, rotate_wilson_frames(grid))
        for cuts in ((nothing, nothing), ([-0.25, 0.25], [-0.25, 0.25]))
            p1, p2 = WilsonWannier2D(frames; wannier_cuts=cuts)
            @test sort(real.(cis.(2π .* p1))) ≈ [-1, 1] atol=1e-12
            @test sort(real.(cis.(2π .* p2))) ≈ [-1, 1] atol=1e-12
        end
    end
    @test_throws r"Wannier sector gap" WilsonWannier2D(grid; wannier_cuts=([0.0, 0.25], nothing))
    @test_throws ArgumentError WilsonWannier2D(grid; wannier_cuts=([0.25, -0.25], nothing))

    # Two identical nonquantized Wannier bands form a rank-two sector. Its
    # nested polarization is the sum, independent of the chosen basis.
    _, broken = statesgrid1D(ssh(0.4, 0.7), N+1, [1]; multimode=:serial)
    ν = Wilson1D(broken[1:N, :, 1])
    for i in 1:N, j in 1:N
        v = kron(broken[i,:,1], broken[j,:,1])
        grid[i,j,1:4,1] = v
        grid[i,j,5:8,2] = v
    end
    for frames in (grid, rotate_wilson_frames(grid))
        p1, p2 = WilsonWannier2D(frames)
        @test cis(2π * only(p1)) ≈ cis(4π*ν) atol=1e-12
        @test cis(2π * only(p2)) ≈ cis(4π*ν) atol=1e-12
    end

    # A rank-two sector separated from a third band must remain invariant
    # even when the occupied gauge mixes states across the sector boundary.
    mixed = zeros(ComplexF64, N, N, 12, 3)
    mixed[:, :, 1:8, 1:2] = grid
    mixed[:, :, 9:12, 3] = conj.(grid[:, :, 1:4, 1])
    for frames in (mixed, rotate_wilson_frames(mixed))
        p1, p2 = WilsonWannier2D(frames; wannier_cuts=([-0.4, 0.0], [-0.4, 0.0]))
        @test cis.(2π .* p1) ≈ [cis(4π*ν), cis(-2π*ν)] atol=1e-12
        @test cis.(2π .* p2) ≈ [cis(4π*ν), cis(-2π*ν)] atol=1e-12
    end

    # Different phases may be grouped deliberately, without defining spurious
    # individual nested bands. The two conjugate copies have total phase zero.
    for i in 1:N, j in 1:N
        grid[i,j,5:8,2] = conj.(grid[i,j,1:4,1])
    end
    p1, p2 = WilsonWannier2D(rotate_wilson_frames(grid); wannier_cuts=([0.0], [0.0]))
    @test p1 ≈ [0] atol=1e-12
    @test p2 ≈ [0] atol=1e-12

    H(k) = ssh(0.4, 0.7)(k)
    for f in (Wilson2D, WilsonWannier2D)
        default = f(H, 11; multimode=:serial)
        explicit = f(H, 11, 11, [1]; multimode=:serial)
        @test first(default) ≈ first(explicit) atol=1e-12
        @test_throws ArgumentError f(H, 1)
    end
end

@testset "Nested Wilson loops: quadrupole symmetry and Chern obstruction" begin
    s0 = Matrix{ComplexF64}(I, 2, 2)
    sx = ComplexF64[0 1; 1 0]; sy = ComplexF64[0 -im; im 0]; sz = ComplexF64[1 0; 0 -1]
    Γ0 = kron(sz, s0); Γ1 = -kron(sy, sx); Γ2 = -kron(sy, sy)
    Γ3 = -kron(sy, sz); Γ4 = kron(sx, s0)
    bbh(t, m) = k -> (t+cos(2π*k[1]))*Γ4 + sin(2π*k[1])*Γ3 +
                      (t+cos(2π*k[2]))*Γ2 + sin(2π*k[2])*Γ1 + m*Γ0
    for (t, expected) in ((0.5, -1), (1.5, 1))
        p1, p2 = WilsonWannier2D(bbh(t, 0.0), 17, 17, 1:2; multimode=:serial)
        @test cis.(2π .* p1) ≈ fill(expected, 2) atol=1e-12
        @test cis.(2π .* p2) ≈ fill(expected, 2) atol=1e-12
    end
    p1, p2 = WilsonWannier2D(bbh(0.5, 0.3), 17, 17, 1:2; multimode=:serial)
    @test all(abs.(imag.(cis.(2π .* p1))) .> 0.1)
    @test all(abs.(imag.(cis.(2π .* p2))) .> 0.1)
    # A Chern band has winding polarization, so no periodic branch exists
    # for the transverse average requested by the nested-polarization API.
    chern(k) = sin(2π*k[1])*sx + sin(2π*k[2])*sy + (-1+cos(2π*k[1])+cos(2π*k[2]))*sz
    @test_throws r"nested polarization winds" WilsonWannier2D(chern, 21, 21, [1]; multimode=:serial)
end

@testset "Nested Wilson loops: sector ranks cannot change" begin
    NX, NY = 9, 17
    frames = zeros(ComplexF64, NX, NY, 4, 2)
    for i in 1:NX, j in 1:NY
        # One phase moves from ~0.1 to ~0.4 across the cut at 0.25,
        # while the other stays near -0.2. No sampled phase sits on the cut.
        a = 0.9 - 0.3sin(π*(i-1)/(NX-1))^2
        frames[i,j,1:2,1] = [sqrt(a), cis(2π*(j-1)/NY)*sqrt(1-a)]
        frames[i,j,3:4,2] = [sqrt(0.8), cis(-2π*(j-1)/NY)*sqrt(0.2)]
    end
    @test_throws r"Wannier sector dimensions change" WilsonWannier2D(frames;
        wannier_cuts=([0.0, 0.25], [0.0]))
    # With no internal boundary this is a valid rank-two sector.
    p1, p2 = WilsonWannier2D(rotate_wilson_frames(frames); wannier_cuts=([0.0], [0.0]))
    @test p1 ≈ [0] atol=1e-12
    @test all(isfinite, p2)
end

@testset "Operator Wilson loops: spin-mixed 2D topological sectors" begin
    sx = ComplexF64[0 1; 1 0]; sy = ComplexF64[0 -im; im 0]; sz = ComplexF64[1 0; 0 -1]
    s0 = Matrix{ComplexF64}(I, 2, 2)
    spin = kron(sz, s0)
    # Two time-reversed Chern blocks with a Rashba term that breaks spin
    # conservation but preserves the energy gap and projected-spin gap.
    function model(k, λ)
        x, y = 2π .* k
        sin(x)*kron(sz,sx) + sin(y)*kron(s0,sy) + (-1+cos(x)+cos(y))*kron(s0,sz) +
            λ*kron(sin(y)*sx-sin(x)*sy, s0)
    end
    for λ in (0.0, 0.3)
        N = 31
        phases = zeros(N, 2)
        gap, egap = Inf, Inf
        for (j,y) in enumerate(range(0, 1; length=N))
            ks, states = statesgrid1Dslice(k -> model(k, λ), 41, 1:2, 1, y; multimode=:serial)
            p, g = sectorwilson(states, spin)
            phases[j,:] = [only(p.negative), only(p.positive)]
            gap = min(gap, g)
            for x in ks
                ev = eigvals(Hermitian(model([x,y], λ)))
                egap = min(egap, ev[3]-ev[2])
            end
            if j == 7
                pr, gr = sectorwilson(rotate_wilson_frames(states), spin)
                @test cis(2π*only(pr.negative)) ≈ cis(2π*only(p.negative)) atol=1e-12
                @test cis(2π*only(pr.positive)) ≈ cis(2π*only(p.positive)) atol=1e-12
                @test gr ≈ g atol=1e-12
            elseif j == 1
                # The Kramers-degenerate Wilson pair still has a well-defined
                # total spin weight, even though individual weights do not.
                result = wilsonweights(states, spin; kgrid=ks)
                @test result.multiplicities == [2]
                @test result.weights ≈ [0] atol=1e-12
            end
        end
        winding = vec(sum(angle.(cis.(2π .* diff(phases; dims=1))); dims=1)) ./ (2π)
        @test abs.(winding) ≈ [1, 1] atol=1e-12
        @test sum(winding) ≈ 0 atol=1e-12
        @test gap > 0.9
        @test egap > 1
    end
end
