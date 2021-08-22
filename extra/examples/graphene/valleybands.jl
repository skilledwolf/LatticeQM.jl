using LinearAlgebra, Plots
# pyplot()

using LatticeQM

println(":: Create lattice...")
lat = Geometries.honeycomb()
println(":: Create hop matrix...")
hops = Operators.graphene(lat; format=:dense)#; mode=:spinhalf)
# Materials.addhaldane!(hops, lat, 0.1; spinhalf=true, mode=:anti)
# addzeeman!(hops, lat, 0.3)

println(":: Create valley operator...")
valley = Operators.valley(lat; spinhalf=false)

println(":: Get bands...")
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, valley)

# Show bands
# save(bands, "playground_graphene/bands.h5")
println(":: Plot bands...")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,240), colorbar=true, markercolor=:PiYG)

mkpath("output")
savefig("output/bands.pdf")
