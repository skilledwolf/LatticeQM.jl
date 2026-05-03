using Test
using LatticeQM
using LatticeQM.Structure: Lattices
using LinearAlgebra: eigvals, norm, dot

# Golden-fixture regression tests for the smallest commensurate twisted bilayer
# graphene cell, N=1 → twist angle ≈ 21.787°. The numbers below were computed
# once with this codebase and are pinned here to detect any silent change in:
#   - the moiré supercell construction (Lattices.twist + foldPC!)
#   - layer assignment via newdimension!("layer", ...)
#   - the inter/intra-layer Slater-Koster hopping in Operators.graphene
#   - the Bloch-sum convention in Hops.fourierphase
#
# All values are deterministic: no RNG, no SCF, just lattice construction +
# diagonalisation of a 28×28 (or 56×56 in spinful) matrix at fixed k.
@testset "TBG: N=1 moiré geometry fixture" begin
    lat = Geometries.honeycomb_twisted(1; fold=true)

    @testset "atom count and dimensions" begin
        # 4 atoms × (1² + 1·1 + 1²) = 4 × 7 sites? Actually for (n,m)=(1,1)
        # the moiré cell holds 4·N_super = 4·7 = 28 atoms (14 per layer).
        @test Lattices.countorbitals(lat) == 28
        @test Lattices.latticedim(lat) == 2
    end

    @testset "layer assignment is 14/14" begin
        @test Lattices.hasdimension(lat, "layer")
        layers = Lattices.extracoordinates(lat, "layer")
        @test sort(unique(layers)) == [0.0, 1.0]
        @test sum(layers .== 0.0) == 14
        @test sum(layers .== 1.0) == 14
    end

    @testset "moiré lattice vectors: |A| = √21, 60° angle preserved" begin
        A = Lattices.getA(lat)
        @test norm(A[:, 1]) ≈ sqrt(21.0) atol=1e-12
        @test norm(A[:, 2]) ≈ sqrt(21.0) atol=1e-12
        # Twist preserves the relative 60° between Bravais vectors.
        cosθ = dot(A[:, 1], A[:, 2]) / (norm(A[:, 1]) * norm(A[:, 2]))
        @test isapprox(cosθ, 0.5; atol=1e-12)
    end

    @testset "twist angle pinned to ≈ 21.787°" begin
        # Precise value: arccos((3·1² + 3·1·1 + 1²/2) / (3·1² + 3·1·1 + 1²))
        #              = arccos(6.5 / 7) ≈ 0.380251... rad ≈ 21.787°.
        @test isapprox(Lattices.twistangle(1), acos(6.5 / 7); atol=1e-14)
        @test isapprox(Lattices.twistangle(1; degrees=true), 21.78678929826181;
                       atol=1e-10)
    end
end

@testset "TBG: N=1 Slater-Koster spectrum fixture" begin
    lat = Geometries.honeycomb_twisted(1; fold=true)
    hops = Operators.graphene(lat; format=:dense, mode=:nospin, tz=0.32)

    @testset "Hops shape and Hermiticity" begin
        @test hopdim(hops) == 28
        # H(k) at the Γ point of the mini-BZ must be Hermitian.
        HG = Matrix(hops([0.0, 0.0]))
        @test size(HG) == (28, 28)
        @test maximum(abs, HG .- HG') < 1e-12
    end

    # Pinned eigenvalues at Γ. tz=0.32 and AB-like stacking break particle-hole
    # symmetry, so the spectrum is not symmetric about zero.
    @testset "eigenvalue snapshot at Γ" begin
        ev = sort(real.(eigvals(Matrix(hops([0.0, 0.0])))))
        # Lowest 4 levels (deep valence).
        @test isapprox(ev[1], -3.5589039818016364; atol=1e-10)
        @test isapprox(ev[2], -2.4429499132468653; atol=1e-10)
        @test isapprox(ev[3], -1.8275362149119176; atol=1e-10)
        @test isapprox(ev[4], -1.7804326240075463; atol=1e-10)
        # Four "near-flat" central levels (n=14,15,16 indices around half-fill).
        @test isapprox(ev[13], -1.0539291875872567; atol=1e-10)
        @test isapprox(ev[14], -1.0005600595025606; atol=1e-10)
        @test isapprox(ev[15],  1.2328013577303936; atol=1e-10)
        @test isapprox(ev[16],  1.2511298250266116; atol=1e-10)
        # Top 4 levels (deep conduction).
        @test isapprox(ev[end-1], 3.000934702452823;  atol=1e-10)
        @test isapprox(ev[end],   3.00099881605841;   atol=1e-10)
    end

    @testset "near-degenerate central pair at K-mini" begin
        # At the corner of the mini Brillouin zone, the 21.787° TBG model has
        # near-degenerate central pairs. Pinning their values catches any shift
        # from changes to interlayer hopping or twist angle.
        ev = sort(real.(eigvals(Matrix(hops([1/3, -1/3])))))
        @test isapprox(ev[14], -0.008387167163543133; atol=1e-10)
        @test isapprox(ev[15],  0.01040073350158507;  atol=1e-10)
        # Two-fold near-degeneracy at the level of 1e-12.
        @test isapprox(ev[14], ev[13]; atol=1e-12)
        @test isapprox(ev[15], ev[16]; atol=1e-12)
    end
end

@testset "TBG: twist invariants" begin
    @testset "twist input lattices are not mutated" begin
        lat = Geometries.honeycomb()
        sp_before = copy(lat.spacecoordinates)
        el_before = copy(lat.extralabels)
        Lattices.twist(lat, lat, 1)
        @test lat.spacecoordinates == sp_before
        @test lat.extralabels      == el_before
    end

    @testset "(n, m) must be coprime" begin
        @test_throws AssertionError Lattices.twist(Geometries.honeycomb(), 2; m=2)
    end

    @testset "moire_supercell returns the legacy √3-rotated cell" begin
        # No `minimal` switch — the function returns the only supported cell.
        # The smaller minimal commensurate moiré cell is geometrically valid
        # but incompatible with the parity-flip-based twist algorithm here.
        S = Lattices.moire_supercell(1, 1)
        @test S == [1 -2; 2 3]
        @test abs(S[1,1]*S[2,2] - S[1,2]*S[2,1]) == 7   # 3n²+3nm+m² for (1,1)
    end
end

# Smoke + Hermiticity for the four-layer twisted constructors. Pre-fix these
# had zero coverage and were calling the buggy in-place code path.
@testset "TBG: four-layer ABBA / ABAB constructors" begin
    for (name, cons) in [("ABBA", Geometries.honeycomb_twisted_ABBA),
                         ("ABAB", Geometries.honeycomb_twisted_ABAB)]
        @testset "$name N=1" begin
            lat = cons(1)
            @test Lattices.latticedim(lat) == 2
            @test Lattices.countorbitals(lat) > 0
            # 4 layers × 2 atoms/cell × cells/layer = countorbitals; cells/layer
            # follows the same legacy 3n²+3nm+m² = 7 for (n=1, m=1).
            @test Lattices.countorbitals(lat) == 4 * 2 * 7
            h = Operators.graphene(lat; format=:dense, mode=:nospin, tz=0.32)
            HK = Matrix(h([0.0, 0.0]))
            @test maximum(abs, HK .- HK') < 1e-10
        end
    end
end

# ============================================================================
# Pre-rewrite fixture set (do not move/edit lightly).
#
# These eigenvalue/lattice-vector snapshots were captured with the in-place
# `twist` recipe (translate, negate, repeat-and-crop, basis-rotate). They are
# pinned here so any rewrite of the construction has to reproduce them
# bit-for-bit (within fp tolerance) before being merged.
#
# Eigenvalues at *fractional* k-points are invariant under any global rotation
# of the basis — so a clean rewrite that constructs the layers in their
# natural orientation (no final basis rotation) must still match these
# spectra exactly.
# ============================================================================

# Helper: spectrum of the graphene tight-binding model (NN + interlayer 0.32t)
# evaluated at fractional k.
function _tbg_spectrum(lat, k)
    h = Operators.graphene(lat; format=:dense, mode=:nospin, tz=0.32)
    sort(real.(eigvals(Matrix(h(k)))))
end

# Strong physical-correctness check: every atom in each layer of a TBG cell
# must have exactly 3 in-plane nearest neighbours at distance ≈ 1 (the
# original honeycomb bond length). Hermiticity + atom counts + bandwidth pass
# even when the construction silently drops bonds, so this NN-coordination
# test is the canonical guard against malformed moiré cells.
function _layer_nn_counts(lat, layer_id; bond=1.0, tol=0.05, shells=2)
    N      = Lattices.countorbitals(lat)
    layers = Lattices.extracoordinates(lat, "layer")
    pos    = Lattices.allpositions(lat)
    A      = Lattices.getA(lat)

    counts = Int[]
    for i in 1:N
        layers[i] == layer_id || continue
        c = 0
        for di in -shells:shells, dj in -shells:shells
            R = A[:, 1] * di + A[:, 2] * dj
            for j in 1:N
                layers[j] == layer_id || continue
                (i == j && di == 0 && dj == 0) && continue
                d = sqrt((pos[1, i] - pos[1, j] - R[1])^2 +
                         (pos[2, i] - pos[2, j] - R[2])^2 +
                         (pos[3, i] - pos[3, j])^2)
                abs(d - bond) < tol && (c += 1)
            end
        end
        push!(counts, c)
    end
    counts
end

@testset "TBG topology: every atom in every layer has 3 in-plane NN" begin
    @testset "honeycomb_twisted N=$N (legacy supercell)" for N in 1:3
        lat = Geometries.honeycomb_twisted(N)
        @test all(==(3), _layer_nn_counts(lat, 0.0))
        @test all(==(3), _layer_nn_counts(lat, 1.0))
    end

    @testset "honeycomb_twisted_ABBA N=1" begin
        lat = Geometries.honeycomb_twisted_ABBA(1)
        for lid in unique(Lattices.extracoordinates(lat, "layer"))
            @test all(==(3), _layer_nn_counts(lat, lid))
        end
    end

    @testset "honeycomb_twisted_ABAB N=1" begin
        lat = Geometries.honeycomb_twisted_ABAB(1)
        for lid in unique(Lattices.extracoordinates(lat, "layer"))
            @test all(==(3), _layer_nn_counts(lat, lid))
        end
    end
end

@testset "TBG fixtures: legacy supercell" begin
    @testset "N=2 (13.17°)" begin
        lat = Geometries.honeycomb_twisted(2)
        @test Lattices.countorbitals(lat) == 76
        @test norm(Lattices.getA(lat)[:, 1]) ≈ sqrt(57.0) atol=1e-12  # |t|² = 3·4 + 6 + 1 = 19; |A| = √(3·19) = √57
        evΓ = _tbg_spectrum(lat, [0.0, 0.0])
        @test isapprox(evΓ[1],  -3.559081532001172; atol=1e-10)
        @test isapprox(evΓ[end], 3.003119080435181; atol=1e-10)
        evK = _tbg_spectrum(lat, [1/3, -1/3])
        @test isapprox(evK[end÷2],    0.0006573224708178428; atol=1e-10)
        @test isapprox(evK[end÷2+1],  0.0030191851800312744; atol=1e-10)
    end

    @testset "N=3 (9.43°)" begin
        lat = Geometries.honeycomb_twisted(3)
        @test Lattices.countorbitals(lat) == 148
        evΓ = _tbg_spectrum(lat, [0.0, 0.0])
        @test isapprox(evΓ[1],  -3.5593781194147773; atol=1e-10)
        @test isapprox(evΓ[end], 3.006649269356615;  atol=1e-10)
    end
end

# Pre-rewrite minimal-supercell fixtures lived here; they are removed because
# `minimal=true` is currently disabled (it produced a malformed rotated layer
# — see the NN-topology test). Once the minimal-cell construction is
# implemented correctly, re-add fixtures pinned against a verified-correct
# baseline (validated by the same NN-coordination invariant).

@testset "TBG fixtures: ABBA / ABAB N=1" begin
    @testset "ABBA N=1" begin
        lat = Geometries.honeycomb_twisted_ABBA(1)
        @test Lattices.countorbitals(lat) == 56
        evΓ = _tbg_spectrum(lat, [0.0, 0.0])
        @test isapprox(evΓ[1],  -3.887921925345088; atol=1e-9)
        @test isapprox(evΓ[end], 3.0576152590759227; atol=1e-9)
    end

    @testset "ABAB N=1" begin
        lat = Geometries.honeycomb_twisted_ABAB(1)
        @test Lattices.countorbitals(lat) == 56
        evΓ = _tbg_spectrum(lat, [0.0, 0.0])
        @test isapprox(evΓ[1],  -3.887918088102941; atol=1e-9)
        @test isapprox(evΓ[end], 3.0576460746843157; atol=1e-9)
    end
end
