using LinearAlgebra, Plots
gr()

using LatticeQM
using LatticeQM.Operators: graphene, gethaldane, valleyoperator, addzeeman!

lat = Geometries2D.honeycomb()
hops = graphene(lat; mode=:spinhalf)
# Materials.addhaldane!(hops, lat, 0.1; spinhalf=true, mode=:anti)
addzeeman!(hops, lat, 0.3)


valley = valleyoperator(lat; spinhalf=true)

ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, valley)

# Show bands
# save(bands, "playground_graphene/bands.h5")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,240), colorbar=true, markercolor=:PiYG)

mkpath("output")
savefig("output/bands.pdf")