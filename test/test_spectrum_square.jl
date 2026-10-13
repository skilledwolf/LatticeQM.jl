using Test
using LatticeQM
using LinearAlgebra: eigvals
import LinearAlgebra
import SparseArrays

# Analytic regression: nearest-neighbour 1-orbital square lattice with
# hopping t = -1 has dispersion
#     ε(k) = -2(cos 2π k₁ + cos 2π k₂)
# in the fractional-k convention used by LatticeQM (Hops.fourierphase uses
# exp(i 2π k·δL)). Bandwidth = 8|t| = 8.
@testset "Spectrum: square NN bands match -2(cos+cos)" begin
    lat = Geometries.square()
    T = Hops()
    Operators.nearestneighbor!(T, lat, -1.0)

    ks = kpath(lat; num_points=200)
    bands = getbands(T, ks)

    analytic = [-2 * (cos(2π * k[1]) + cos(2π * k[2])) for k in eachcol(ks.points)]
    @test maximum(abs, bands.bands[1, :] .- analytic) < 1e-12
end

# Sanity: bandwidth and band edges over a uniform BZ grid.
@testset "Spectrum: square bandwidth and edges" begin
    lat = Geometries.square()
    T = Hops()
    Operators.nearestneighbor!(T, lat, -1.0)

    # Sweep a small uniform grid; min and max sit at Γ and M.
    # k-vectors must match latticedim(square)=2; passing a 3-vector errors
    # in the Bloch dot-product.
    ks = [[kx, ky] for kx in range(0, 1; length=21), ky in range(0, 1; length=21)]
    energies = [real(only(eigvals(Matrix(T(k))))) for k in ks]
    @test isapprox(minimum(energies), -4.0; atol=1e-12)
    @test isapprox(maximum(energies),  4.0; atol=1e-12)
end

@testset "Spectrum: chemical potential honors filling and k weights" begin
    bands = reshape([-3.0, -1.0, 1.0, 3.0], 1, :)
    ks = zeros(2, 4)
    mesh = Structure.Mesh(ks, [0.1, 0.15, 0.25, 0.5])
    path = Structure.Paths.DiscretePath(ks, [1.0, 4.0], ["start", "end"], collect(1.0:4.0))
    for (grid, targets) in ((ks, [0.0, 0.25, 0.5, 0.75, 1.0]),
                            (path, [0.0, 0.25, 0.5, 0.75, 1.0]),
                            (mesh, [0.0, 0.1, 0.25, 0.5, 1.0]))
        for target in targets
            μ = Spectrum.chemicalpotential(bands, grid, target)
            @test Spectrum.getfilling(bands, grid; μ) ≈ target
        end
        for T in (0.05, 2.0), target in (0.0, 1e-6, 0.02, 0.8, 1-1e-6, 1.0)
            μ = Spectrum.chemicalpotential(bands, grid, target; T)
            @test Spectrum.getfilling(bands, grid; μ, T) ≈ target atol=1e-12
        end
    end
    flatbands = reshape(repeat([-1.0, 1.0]; inner=5000), 1, :)
    uniformmesh = Structure.Mesh(zeros(2, length(flatbands)))
    μ = Spectrum.chemicalpotential(flatbands, uniformmesh, 0.5)
    @test Spectrum.getfilling(flatbands, uniformmesh; μ) ≈ 0.5
    @test_throws ArgumentError Spectrum.chemicalpotential(bands, ks, -0.1)
    @test_throws ArgumentError Spectrum.chemicalpotential(bands, ks, 0.5; T=-1)
end

@testset "Spectrum: multi-frequency LDOS is independent of k partitions" begin
    energies = [-1.0, 0.7]
    H = k -> Matrix(LinearAlgebra.Diagonal(energies))
    ωs = [-0.5, 0.0, 0.3]
    Γ = 0.2
    expected = [sum(Γ / (π * ((ω-e)^2 + Γ^2)) for ω in ωs) / length(ωs) for e in energies]
    for nk in (1, 2, 8), multimode in (:serial, :threaded)
        got = Spectrum.ldos(H, zeros(2, nk), ωs; Γ, multimode, format=:dense, hidebar=true)
        @test got ≈ expected
    end
end

@testset "Spectrum: sparse density forwards eigensolver keywords" begin
    H = k -> SparseArrays.spdiagm(0 => [-4.0, -3.0, -2.0, -1.0, 1.0, 2.0])
    ks = zeros(2, 3)
    expected = Spectrum.density(H, ks; format=:dense, multimode=:serial, hidebar=true)
    actual = Spectrum.density(H, ks; format=:sparse, num_bands=4, which=:SR,
                              multimode=:serial, hidebar=true)
    @test actual ≈ expected
end

@testset "Spectrum: Fermi surface lattice coordinates" begin
    lat = Geometries.square()
    H = k -> reshape([cos(2π*k[1]) + cos(2π*k[2])], 1, 1)
    ks, density = Spectrum.fermisurfacedensity(H, 0.0; lat, num_points=3,
        broadening=0.2, multimode=:serial, hidebar=true)
    @test size(ks) == (3, 9)
    @test size(density) == (1, 9)
    @test all(isfinite, ks)
    @test all(isfinite, density)
end
