using LatticeQM
using LatticeQM.Operators: graphene

println("Generating lattice geometry...")
@time lat = Geometries2D.honeycomb_twisted(30)

println("Generating hops...")
@time  (hops = graphene(lat; format=:sparse, mode=:nospin, tz=0.12))

println("Done.")
