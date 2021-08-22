using LinearAlgebra, Plots

using LatticeQM
import LatticeQM.Structure

lat = Geometries.honeycomb()
H = Operators.graphene(lat; mode=:spinhalf)
H = DenseHops(H)

energies, dos = Spectrum.getdos(H, -3.1, 3.1; klin=800, format=:dense, Γ=0.005)

plot(energies, dos)

mkpath("output")
savefig("output/dos.pdf")