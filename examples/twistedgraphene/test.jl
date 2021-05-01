using LatticeQM

println("Generating lattice geometry...")
@time lat = Structure.Geometries.honeycomb_twisted(30)

println("Generating hops...")
@time  (hops = Operators.graphene(lat; format=:sparse, mode=:nospin, tz=0.12))

println("Done.")
