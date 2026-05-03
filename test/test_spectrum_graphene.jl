using Test
using LatticeQM

# Fixture: analytic graphene conduction band ε_+(k)
include("analyticbands.jl")
using .analyticbands: grapheneconductionband

# Regression: nearest-neighbour graphene tight-binding bands must reproduce the
# analytic two-band result ε_±(k) = ±|t|·|f(k)| to numerical precision along
# the high-symmetry k-path. This pins down the joint behaviour of:
#   - Geometries.honeycomb()
#   - Operators.nearestneighbor!
#   - Operators.setfilling!
#   - kpath / Spectrum.getbands
@testset "Spectrum: graphene NN bands match analytic formula" begin
    lat = Geometries.honeycomb()
    T = Hops()
    Operators.nearestneighbor!(T, lat)
    Operators.setfilling!(T, 0.5; nk=100)

    ks = kpath(lat; num_points=200)
    bands = getbands(T, ks)
    upper_analytic = map(grapheneconductionband, eachcol(ks.points))
    lower_analytic = -upper_analytic

    mae_lower = sum(abs, bands.bands[1, :] .- lower_analytic)
    mae_upper = sum(abs, bands.bands[2, :] .- upper_analytic)

    @test mae_lower + mae_upper < 1e-12
end
