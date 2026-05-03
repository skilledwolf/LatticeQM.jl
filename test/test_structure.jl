using Test
using LatticeQM
using LatticeQM.Structure: Lattices
using LinearAlgebra: I, norm, det

# Sanity checks for Lattice constructors and dimensional accessors. These are
# fast and catch any change to the Lattice schema (basis layout, default
# extracoordinates, latticedim, etc.).
@testset "Structure: 2D Geometries shapes" begin
    @testset "honeycomb()" begin
        lat = Geometries.honeycomb()
        @test Lattices.spacedim(lat) == 3       # Lattice always embeds in ≥ d
        @test Lattices.latticedim(lat) == 2
        @test Lattices.countorbitals(lat) == 2
        @test Lattices.hasdimension(lat, "sublattice")
        @test Lattices.extracoordinates(lat, "sublattice") == [0.0, 1.0]
    end

    @testset "square()" begin
        lat = Geometries.square()
        @test Lattices.latticedim(lat) == 2
        @test Lattices.countorbitals(lat) == 1
        # Square Bravais vectors are orthogonal of length 1.
        A = Lattices.getA(lat)
        @test norm(A[:, 1]) ≈ 1.0
        @test norm(A[:, 2]) ≈ 1.0
        @test isapprox(A[:, 1]' * A[:, 2], 0; atol=1e-12)
    end

    @testset "triangular()" begin
        lat = Geometries.triangular()
        @test Lattices.latticedim(lat) == 2
        @test Lattices.countorbitals(lat) == 1
        # Triangular Bravais vectors of length 1, internal angle 60°.
        A = Lattices.getA(lat)
        @test norm(A[:, 1]) ≈ 1.0
        @test norm(A[:, 2]) ≈ 1.0
        cosθ = (A[:, 1]' * A[:, 2]) / (norm(A[:, 1]) * norm(A[:, 2]))
        @test isapprox(cosθ, 0.5; atol=1e-12)   # 60°
    end

    @testset "triangular_supercell() carries sublattice labels 1,2,3" begin
        lat = Geometries.triangular_supercell()
        @test Lattices.countorbitals(lat) == 3
        @test sort(Lattices.extracoordinates(lat, "sublattice")) == [1.0, 2.0, 3.0]
    end

    @testset "honeycomb_bilayer / AA / AB / BA constructors" begin
        for cons in (Geometries.honeycomb_bilayer,
                     Geometries.honeycomb_AA,
                     Geometries.honeycomb_AB,
                     Geometries.honeycomb_BA)
            lat = cons()
            @test Lattices.countorbitals(lat) == 4
            @test Lattices.hasdimension(lat, "sublattice")
        end
    end
end

# The reciprocal lattice satisfies A^T B = I (per LatticeQM's convention,
# B = A (A^T A)^{-1}). Catches any change in the dual-lattice convention.
@testset "Structure: dual basis satisfies A^T B = I" begin
    for lat in (Geometries.honeycomb(),
                Geometries.square(),
                Geometries.triangular(),
                Geometries.triangular_supercell())
        A = Lattices.getA(lat)
        B = Lattices.getB(lat)
        d = Lattices.latticedim(lat)
        @test A[:, 1:d]' * B[:, 1:d] ≈ Matrix{Float64}(I, d, d)
    end
end

# kpath produces a path with the expected total number of sample points and the
# correct k-vector dimension.
@testset "Structure: kpath shape" begin
    lat = Geometries.honeycomb()
    ks = kpath(lat; num_points=200)
    @test size(ks)[2] == 200
    @test size(ks.points, 1) == 2   # 2D BZ
end
