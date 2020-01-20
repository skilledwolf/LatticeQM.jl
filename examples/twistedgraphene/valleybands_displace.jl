using LinearAlgebra, Plots
gr()

using LatticeQM
using LatticeQM.Operators: graphene, valleyoperator
using LatticeQM.Operators: addinterlayerbias!, setfilling!, gethaldane, addrashba!, addzeeman!, valleyoperator

println("Generating lattice geometry...")
lat = Geometries2D.honeycomb_twisted(11)
Geometries2D.smoothdisplaceZ!(lat, 0.045)

println("Generating lattice operators...")
valley = valleyoperator(lat; spinhalf=false)
hops = graphene(lat; format=:sparse, mode=:nospin, tz=0.32)
setfilling!(hops, lat, 0.5+0.0/hopdim(hops); nk=9)
# addinterlayerbias!(hops, lat, 0.05)
# addhaldane!(hops, lat, 0.1; spinhalf=true, mode=:anti)
# addzeeman!(hops, lat, 0.05)

println("Calculating bands...")
ks = kpath(lat; num_points=100)
bands = getbands(hops, ks, valley; format=:sparse, num_bands=36)

# Show bands
# save(bands, "playground_graphene/bands.h5")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,240), colorbar=true, markercolor=:PiYG)

mkpath("output"); savefig("output/valleybands_displace.pdf")