using Test
using LatticeQM
using LinearAlgebra: eigvals

# Analytic regression: nearest-neighbour 1-orbital triangular lattice with
# hopping t = -1 has six neighbours at fractional vectors ±(1,0), ±(0,1),
# ±(1,-1), giving
#     ε(k) = -2[cos 2π k₁ + cos 2π k₂ + cos 2π(k₁ - k₂)].
# Band edges: ε_min = -6 at Γ; ε_max = +3 at K (k = (1/3, -1/3) and equivalents).
@testset "Spectrum: triangular NN bands match analytic" begin
    lat = Geometries.triangular()
    T = Hops()
    Operators.nearestneighbor!(T, lat, -1.0)

    ks = kpath(lat; num_points=200)
    bands = getbands(T, ks)

    analytic = [-2 * (cos(2π * k[1]) +
                      cos(2π * k[2]) +
                      cos(2π * (k[1] - k[2])))
                for k in eachcol(ks.points)]
    @test maximum(abs, bands.bands[1, :] .- analytic) < 1e-12
end

@testset "Spectrum: triangular band edges (Γ→-6, K→+3)" begin
    lat = Geometries.triangular()
    T = Hops()
    Operators.nearestneighbor!(T, lat, -1.0)

    # k must be 2-component to match latticedim(triangular)=2.
    @test real(only(eigvals(Matrix(T([0.0, 0.0]))))) ≈ -6.0
    @test real(only(eigvals(Matrix(T([1/3, -1/3]))))) ≈ 3.0
end
