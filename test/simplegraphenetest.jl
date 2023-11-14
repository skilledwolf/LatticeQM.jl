module simplegraphenetest
using Test
using LatticeQM

include("analyticbands.jl")
using .analyticbands: grapheneconductionband

function graphenebandMAE()
    lat = Geometries.honeycomb()
    T = Hops()
    Operators.nearestneighbor!(T, lat)
    Operators.setfilling!(T, 0.5; nk=100)

    ks = kpath(lat; num_points=200)
    bands = getbands(T, ks)
    conductionband = map(grapheneconductionband, eachcol(ks.points))
    valenceband = - conductionband

    maeconduction = sum(broadcast(abs, bands.bands[1, :] - valenceband))
    maevalence = sum(broadcast(abs, bands.bands[2, :] - conductionband))
    mae = maeconduction + maevalence
    return mae
end

mae = graphenebandMAE()
@test (mae < 1e-12)
end