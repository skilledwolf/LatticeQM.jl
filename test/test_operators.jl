using Test
using LatticeQM
using LinearAlgebra: I, tr, Diagonal

# Spin operators on a 2-orbital honeycomb sit in a (2N)×(2N) = 4×4 space.
# They must satisfy the spin-1/2 Pauli algebra block-by-block.
@testset "Operators: spin operators on honeycomb" begin
    lat = Geometries.honeycomb()
    N = 2
    sx = getoperator(lat, "sx")
    sy = getoperator(lat, "sy")
    sz = getoperator(lat, "sz")

    @testset "shape" begin
        @test size(sx) == (2N, 2N)
        @test size(sy) == (2N, 2N)
        @test size(sz) == (2N, 2N)
    end

    @testset "Hermiticity" begin
        @test sx ≈ sx'
        @test sy ≈ sy'
        @test sz ≈ sz'
    end

    @testset "involution: sᵢ² = I (in the 2N×2N space)" begin
        Imat = Matrix{ComplexF64}(I, 2N, 2N)
        @test sx * sx ≈ Imat
        @test sy * sy ≈ Imat
        @test sz * sz ≈ Imat
    end

    @testset "commutation: [sx, sy] = 2i sz (and cyclic)" begin
        @test sx * sy - sy * sx ≈ 2im * sz
        @test sy * sz - sz * sy ≈ 2im * sx
        @test sz * sx - sx * sz ≈ 2im * sy
    end

    @testset "name-aliasing: 'SX' / 'sx' / 'MX' all return same operator" begin
        @test getoperator(lat, "SX") ≈ sx
        @test getoperator(lat, "MX") ≈ sx
    end

    @testset "spin up/down projectors sum to identity and are orthogonal" begin
        sup = getoperator(lat, "spinup")
        sdn = getoperator(lat, "spindown")
        Imat = Matrix{ComplexF64}(I, 2N, 2N)
        @test sup + sdn ≈ Imat
        @test sup * sdn ≈ zeros(2N, 2N)
        @test sup * sup ≈ sup       # idempotent
    end
end

# Sublattice projectors must be orthogonal idempotents that partition the
# Hilbert space.
@testset "Operators: sublattice projectors on honeycomb" begin
    lat = Geometries.honeycomb()
    PA = getoperator(lat, "sublatticeA")
    PB = getoperator(lat, "sublatticeB")
    @test size(PA) == (2, 2)        # honeycomb has N=2 orbitals; spinless here
    @test PA + PB ≈ Matrix{Float64}(I, 2, 2)
    @test PA * PB ≈ zeros(2, 2)
    @test PA * PA ≈ PA
    @test tr(PA) ≈ 1
    @test tr(PB) ≈ 1
end

# Sublattice-spin projectors live in (2N)×(2N) and project onto e.g. all
# A-sublattice spin states.
@testset "Operators: sublattice-spin projectors" begin
    lat = Geometries.honeycomb()
    PAs = getoperator(lat, "sublatticeAspin")
    PBs = getoperator(lat, "sublatticeBspin")
    @test size(PAs) == (4, 4)
    @test PAs + PBs ≈ Matrix{Float64}(I, 4, 4)
    @test PAs * PBs ≈ zeros(4, 4)
    @test tr(PAs) ≈ 2     # 1 sublattice site × 2 spin states
end

# trace and expval on Hops; uses the onsite (zero-key) sector.
@testset "Operators: trace and expval on identity Hops" begin
    Imat = Matrix{ComplexF64}(I, 4, 4)
    H = Hops([0, 0] => Imat)
    @test Operators.trace(H) ≈ 4
    @test Operators.expval(H) ≈ tr(Imat)   # sum of all entries — onsite only here
end

# Regression: Spectrum/Hops current operators J_α(k) = ∂H(k)/∂k_α must be
# Hermitian for any Hermitian H. The Hops-input branch of getcurrentoperators
# is the documented "preferred path" for tight-binding Hamiltonians, so this
# is the right one to pin. Catches sign-flip / position-convention drift in
# `currentoperators!` and the dense/sparse helpers it dispatches to.
@testset "Operators: getcurrentoperators(Hops) is Hermitian" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    Js = Operators.getcurrentoperators(lat, H)

    # One operator per lattice direction (graphene = 2). The function
    # previously returned `spacedim(lat) = 3` operators with the third one
    # identically zero; current.jl now returns only the meaningful ones.
    @test length(Js) == LatticeQM.Structure.Lattices.latticedim(lat)

    for k in ([0.0, 0.0], [1/3, -1/3], [0.5, 0.5], [-0.2, 0.7])
        for J in Js
            Jk = Matrix(J(k))
            @test maximum(abs, Jk .- Jk') < 1e-10
        end
    end
end
