using Test
using LatticeQM
using LinearAlgebra: eigvals

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
