using LinearAlgebra
using Plots

using LatticeQM


############################################################
# Build geometries
############################################################
N=10
lat_armchair = Structure.Lattices.reduceto1D(Geometries2D.honeycomb(), [[1, 1] [N, -N]])
lat_zigzag   = Structure.Lattices.reduceto1D(Geometries2D.honeycomb(), [[1, 0] [0, N]])

p1 = plot(lat_armchair; supercell=[0:20])
p2 = plot(lat_zigzag; supercell=[0:15])
p = plot(p1,p2, layout=(2,1))
savefig(p, "output/ribbons_geometry.pdf")

############################################################
# Bands (armchair)
############################################################
hops = Operators.graphene(lat_armchair)
ks = kpath(lat_armchair; num_points=800);
bands1 = getbands(hops, ks)
p1 = plot(bands1, 0; ylabel="\$\\varepsilon/t\$", size=(330,260))

############################################################
# Bands (zigzag)
############################################################
hops = Operators.graphene(lat_zigzag)
ks = kpath(lat_zigzag; num_points=800);
bands2 = getbands(hops, ks)
p2 = plot(bands2, 0; ylabel="\$\\varepsilon/t\$", size=(330,260))

p = plot(p1,p2, layout=(2,1))
savefig(p, "output/ribbons_bands.pdf")