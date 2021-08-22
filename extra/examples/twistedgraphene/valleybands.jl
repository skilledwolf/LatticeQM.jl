using LinearAlgebra, Plots
gr()

using LatticeQM

println("Generating lattice geometry...")
@time lat = Geometries.honeycomb_twisted(5)

println("Generating valley operator...")
@time valley = Operators.valley(lat; spinhalf=false)

println("Generating Hamiltonian operator...")
@time hops = Operators.graphene(lat; format=:sparse, mode=:nospin, tz=0.32)

println("Setting filling...")
@time Operators.setfilling!(hops, 0.5+0.0/hopdim(hops); nk=9)
# addinterlayerbias!(hops, lat, 0.05)
# addhaldane!(hops, lat, 0.1; spinhalf=true, mode=:anti)
# addzeeman!(hops, lat, 0.05)

println("Calculating bands...")
ks = kpath(lat; num_points=100)
@time bands = getbands(hops, ks, valley; format=:sparse, num_bands=36)

# Show bands
# save(bands, "playground_graphene/bands.h5")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,240), colorbar=true, markercolor=:PiYG)

mkpath("output"); savefig("output/valleybands.pdf")
