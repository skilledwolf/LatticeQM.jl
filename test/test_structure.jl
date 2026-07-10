using Test
using LatticeQM
using LatticeQM.Structure: Lattices
using LinearAlgebra: I, norm, det, dot

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

# ---------------------------------------------------------------------------
# Wigner–Seitz folding and neighbor-cell search on non-hexagonal bases.
# The old implementations built their candidate sets from the raw basis
# (single shortest shell / fixed integer box), so unequal-length or sheared
# bases silently missed folds and true nearest cells.
# ---------------------------------------------------------------------------
@testset "Structure: foldcell! on unequal-length and sheared bases" begin
    # rectangular diag(1,2): the long axis used to be left completely unfolded
    pts = Float64[0.2 0.0 -0.7; 5.3 -3.49 1.0]
    Lattices.foldcell!(pts, [1.0 0.0; 0.0 2.0])
    A = [1.0 0.0; 0.0 2.0]
    M = A' * A
    for p in eachcol(pts)
        # WS condition: p is closer to the origin than to any lattice point
        for n1 in -2:2, n2 in -2:2
            (n1 == 0 && n2 == 0) && continue
            g = [n1, n2]
            @test dot(p, M, p) <= dot(p - g, M, p - g) + 1e-9
        end
    end

    # sheared basis [1 6; 0 1] (unimodular shear of the square lattice):
    # folding must reproduce square-lattice behavior
    pts2 = Float64[7.3; 0.9;;]
    Lattices.foldcell!(pts2, [1.0 6.0; 0.0 1.0])
    cart = [1.0 6.0; 0.0 1.0] * pts2[:, 1]
    @test maximum(abs, cart) <= 0.5 + 1e-9   # square-lattice WS cell
end

@testset "Structure: getneighborcells finds true nearest cells (sheared basis)" begin
    A = [1.0 6.0; 0.0 1.0]   # same lattice as the square lattice
    cells = Lattices.getneighborcells(A, 1; halfspace=false, innerpoints=false,
                                      excludeorigin=true)
    dists = sort([norm(A * v) for v in cells])
    @test length(cells) == 4                  # square lattice: 4 nearest cells
    @test all(isapprox.(dists, 1.0; atol=1e-9))
    @test [-6, 1] in cells || [6, -1] in cells  # the shell the fixed box missed

    # halfspace dedup keeps exactly one of each ± pair
    half = Lattices.getneighborcells(A, 1; halfspace=true, innerpoints=false,
                                     excludeorigin=true)
    @test length(half) == 2
    for v in half
        @test !(-v in half)
    end

    # regression guard: reduced-basis search must not change hexagonal results
    lat = Geometries.honeycomb()
    c1 = Lattices.getneighborcells(lat, 1; halfspace=false, innerpoints=true,
                                   excludeorigin=false)
    @test length(c1) == 7                     # origin + 6 hexagonal neighbors
end

# ---------------------------------------------------------------------------
# Label bookkeeping: extra-coordinate labels must stay paired with their rows,
# and derived lattices must not alias the input's mutable containers.
# ---------------------------------------------------------------------------
@testset "Structure: cropcircle preserves label→row pairing" begin
    # Two extra dimensions whose Dict hash order differs from row order is
    # exactly the case that used to come back swapped.
    lat = Geometries.honeycomb()                       # has "sublattice" (row 1)
    Lattices.newdimension!(lat, "zheight", fill(7.0, (1, Lattices.countorbitals(lat))))

    circ = Lattices.circleregion(lat, 3.0)
    @test Lattices.countorbitals(circ) > 0
    subl = Lattices.extracoordinates(circ, "sublattice")
    zz = Lattices.extracoordinates(circ, "zheight")
    @test sort(unique(subl)) ⊆ [0.0, 1.0]              # sublattice values, not 7.0
    @test all(zz .== 7.0)                              # zheight values, not 0/1
end

@testset "Structure: superlattice does not alias the input lattice" begin
    lat = Geometries.honeycomb()
    n_extra_before = length(lat.extralabels)
    sup = Lattices.superlattice(lat, [2, 2])
    Lattices.newdimension!(sup, "layer", zeros(1, Lattices.countorbitals(sup)))
    # the input lattice must be untouched (used to gain a phantom "layer" key)
    @test length(lat.extralabels) == n_extra_before
    @test !Lattices.hasdimension(lat, "layer")
    @test Lattices.hasdimension(sup, "layer")
    # and mutating the supercell's specialpoints must not touch the input's
    @test !(sup.specialpoints === lat.specialpoints)
end

@testset "Structure: superlattice maps special k-points to the supercell BZ" begin
    lat = Geometries.honeycomb()
    S = [2 0; 0 2]
    sup = Lattices.superlattice(lat, S)
    B = Lattices.getB(lat)
    Bsup = Lattices.getB(sup)
    d = Lattices.latticedim(lat)
    for name in keys(lat.specialpoints.points)
        p = lat.specialpoints.points[name]
        length(p) == d || continue
        # same Cartesian k-point in both descriptions
        @test B[:, 1:d] * p ≈ Bsup[:, 1:d] * sup.specialpoints.points[name] atol=1e-12
    end
end

@testset "Structure: fillregion builds a labelled 0D flake" begin
    lat = Geometries.honeycomb()
    flake = Lattices.fillregion(lat, p -> norm(p) < 2.5)   # used to throw UndefVarError
    @test Lattices.latticedim(flake) == 0
    @test Lattices.countorbitals(flake) > 6
    @test Lattices.hasdimension(flake, "sublattice")
    @test sort(unique(Lattices.extracoordinates(flake, "sublattice"))) ⊆ [0.0, 1.0]
end

# ---------------------------------------------------------------------------
# Previously throw-on-first-call utilities and geometry helpers.
# ---------------------------------------------------------------------------
@testset "Structure: copy/mirrorZ/mergelattices/rotation3D/signedangle" begin
    lat = Geometries.honeycomb_AB()

    # copy is deep enough: mutating the copy leaves the original alone
    latc = copy(lat)
    latc.spacecoordinates[3, 1] += 1.0
    @test lat.spacecoordinates[3, 1] != latc.spacecoordinates[3, 1]

    # mirrorZ (non-mutating) used to throw MethodError (no copy(::Lattice))
    lam = Lattices.mirrorZ(lat)
    @test lam.spacecoordinates[3, :] ≈ -lat.spacecoordinates[3, :]
    @test lam.spacecoordinates !== lat.spacecoordinates

    # mergelattices (non-mutating) used to throw MethodError
    n = Lattices.countorbitals(lat)
    merged = Lattices.mergelattices(lat, lat)
    @test Lattices.countorbitals(merged) == 2n
    @test Lattices.countorbitals(lat) == n            # input untouched

    # rotation3D must be a proper rotation (old version was symmetric ⇒ not one)
    for (θ, ax) in ((0.3, [0, 0, 1.0]), (1.1, [1.0, 2.0, -0.5]), (-0.7, [0, 1.0, 0]))
        R = LatticeQM.Structure.rotation3D(θ, ax)
        @test R * R' ≈ Matrix(1.0I, 3, 3) atol=1e-12
        @test det(R) ≈ 1.0 atol=1e-12
        @test R * ax ./ norm(ax) ≈ ax ./ norm(ax) atol=1e-12   # axis is fixed
    end
    @test LatticeQM.Structure.rotation3D(0.3) * [1, 0, 0] ≈ [cos(0.3), sin(0.3), 0.0] atol=1e-14

    # signedangle used to reference an undefined variable
    @test LatticeQM.Structure.signedangle([1.0, 0, 0], [0.0, 1, 0]) ≈ π / 2
    @test LatticeQM.Structure.signedangle([0.0, 1, 0], [1.0, 0, 0]) ≈ -π / 2
end

@testset "Structure: getpath handles tick-coincident samples" begin
    # 3 collinear equidistant points; num=5 lands a sample exactly on the
    # middle tick, which used to be pushed twice (6 columns for num=5).
    pts = [0.0 1.0 2.0; 0.0 0.0 0.0]
    path, ticks, s = LatticeQM.Structure.Paths.getpath(pts; num=5)
    @test size(path, 2) == 5
    @test length(s) == 5
    @test ticks ≈ [0.0, 1.0, 2.0]
end

@testset "Structure: regulargrid returns the full requested grid" begin
    @test size(LatticeQM.Structure.regulargrid(nk=1000, dim=3), 2) == 1000
    @test size(LatticeQM.Structure.regulargrid(nk=100, dim=2), 2) == 100
end
