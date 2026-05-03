using Test
using LatticeQM
using LatticeQM.TightBinding: hermitianize!, ishermitian
using SparseArrays: issparse

# Type-system / constructor coverage for Hops and its dense/sparse variants.
@testset "TightBinding: Hops constructors and types" begin
    @testset "empty Hops + addition by setindex!" begin
        H = Hops()
        H[[0, 0]] = ComplexF64[1.0 0.0; 0.0 -1.0]
        @test length(H) == 1
        @test hopdim(H) == 2
        @test H[[0, 0]] == ComplexF64[1.0 0.0; 0.0 -1.0]
    end

    @testset "Hops(::Pair) constructor" begin
        H = Hops([0, 0] => ComplexF64[0 1; 1 0],
                 [1, 0] => ComplexF64[0 0; 1 0])
        @test length(H) == 2
        @test hopdim(H) == 2
    end

    @testset "DenseHops vs SparseHops dispatch" begin
        M = ComplexF64[0 1; 1 0]
        Hd = DenseHops([0, 0] => M)
        Hs = SparseHops([0, 0] => M)
        @test Hd isa DenseHops
        @test Hs isa SparseHops
        @test issparse(Hs)
        @test !issparse(Hd)
        # dense ↔ sparse round-trip preserves contents
        @test dense(Hs)[[0, 0]] ≈ Hd[[0, 0]]
        @test sparse(Hd)[[0, 0]] ≈ Hs[[0, 0]]
    end
end

# Algebraic structure of the Hops container.
@testset "TightBinding: arithmetic on Hops" begin
    H1 = Hops([0, 0] => ComplexF64[1 0; 0 1])
    H2 = Hops([0, 0] => ComplexF64[0 1; 1 0])

    @testset "addition is commutative and componentwise" begin
        @test (H1 + H2)[[0, 0]] ≈ H1[[0, 0]] + H2[[0, 0]]
        @test (H2 + H1)[[0, 0]] ≈ (H1 + H2)[[0, 0]]
    end

    @testset "scalar multiplication" begin
        @test (2 * H1)[[0, 0]] ≈ 2 .* H1[[0, 0]]
        @test (H1 * 2)[[0, 0]] ≈ 2 .* H1[[0, 0]]
    end

    @testset "addhops!(H, H') mutates and returns H" begin
        H = Hops([0, 0] => ComplexF64[1 0; 0 1])
        out = addhops!(H, Hops([0, 0] => ComplexF64[0 0; 0 1]))
        @test out === H
        @test H[[0, 0]] ≈ ComplexF64[1 0; 0 2]
    end
end

# Hermiticity helpers: must consider both H[R] and H[-R]† blocks.
@testset "TightBinding: hermitianize!/ishermitian" begin
    # Onsite-only block must already be Hermitian.
    Hgood = Hops([0, 0] => ComplexF64[1 1im; -1im 2])
    @test ishermitian(Hgood)

    # Construct a non-Hermitian H, hermitianize, then check.
    Hbad = Hops([0, 0]  => ComplexF64[0 1; 0 0],
                [1, 0]  => ComplexF64[0 0; 0 1],
                [-1, 0] => ComplexF64[0 0; 0 0])
    @test !ishermitian(Hbad)
    hermitianize!(Hbad)
    @test ishermitian(Hbad)
end

# Sanity: nearestneighbor! on honeycomb produces the expected number of
# hopping vectors (one per A↔B nearest-neighbour direction). Higher-level
# correctness is covered by the analytic-spectrum regression tests in
# test_spectrum_graphene.jl etc.
@testset "TightBinding: nearestneighbor! shape on honeycomb" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)

    @test hopdim(H) == 2
    @test length(H) >= 1
    # The onsite block (R=0) carries the A↔B hopping (within unit cell).
    @test haskey(H, [0, 0])
end

# Regression: the hopping rule used by nearestneighbor! must be a *concrete
# struct* with the same type across all calls (parameterised only by typeof(t0)),
# so that the downstream gethops/hops!/hoppingmatrix! pipeline compiles once
# per t0-type and the precompile workload in src/LatticeQM.jl can cache it.
# A previous closure-based implementation produced a fresh anonymous-function
# type on every call, defeating precompilation.
@testset "TightBinding: nearestneighbor! hopping rule is a stable concrete type" begin
    lat = Geometries.honeycomb()
    H1, H2 = Hops(), Hops()
    Operators.nearestneighbor!(H1, lat, -1.0)
    Operators.nearestneighbor!(H2, lat, -2.5)
    # Identical hopping matrices up to the t0 scaling (catches any drift
    # in the struct → addhops! plumbing).
    @test H1[[0, 0]] * 2.5 ≈ H2[[0, 0]]

    # The rule itself is constructible and reusable.
    rule = Operators.DistanceWindowHopping(1.0, 0.01, -1.0)
    @test rule isa Operators.DistanceWindowHopping{Float64}
    @test rule(1.0)   == -1.0    # in-window
    @test rule(0.5)   == 0.0     # below window
    @test rule(2.0)   == 0.0     # above window

    # The struct is parametric over t0's type and stays type-stable: in
    # particular zero(T) (not Float64 0.0) is returned for the
    # out-of-window branch. Note: passing complex t0 to nearestneighbor!
    # would produce a non-Hermitian Hamiltonian downstream — the struct
    # itself is correct, but the gethops convention is intended for
    # symmetric real-valued hopping rules. Complex hopping should be
    # introduced via Peierls-type routines instead.
    crule = Operators.DistanceWindowHopping(1.0, 0.01, -1.0 + 0.5im)
    @test crule isa Operators.DistanceWindowHopping{ComplexF64}
    @test crule(1.0) === -1.0 + 0.5im
    @test crule(0.5) === 0.0 + 0.0im
end
