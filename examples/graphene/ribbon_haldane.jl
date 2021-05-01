using LinearAlgebra
using Plots

using LatticeQM


############################################################
# Build geometry
############################################################
N=30
lat1 = Structure.Lattices.reduceto1D(Structure.Geometries.honeycomb(), [[1, 0] [0, N]])
ks = kpath(lat1; num_points=500)

pos = Operators.positionalong(lat1, Structure.basis(lat1,2))


############################################################
# Bands
############################################################
hops1 = Operators.graphene(lat1; cellrange=1, mode=:nospin, tz=0.3, ℓinter=0.08, ℓintra=0.08)
Operators.addhaldane!(hops1, lat1, 0.3)

bands1 = getbands(hops1, ks, pos)
p = plot(bands1, 1; ylabel="\$\\varepsilon/t\$", colorbar_title="transverse position", markersize=1.25, size=(500,300), cquantile=1.0, colorbar=true, csymmetric=false)
savefig(p, "output/ribbon_haldane_bands.pdf")