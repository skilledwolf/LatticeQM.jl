using LinearAlgebra, Plots
gr()

using LatticeQM
using LatticeQM.Operators: graphene, gethaldane, valleyoperator, addzeeman!

lat = Geometries2D.honeycomb()
hops = graphene(lat; mode=:spinhalf)

ks = LatticeQM.Utils.regulargrid(nk=400^2)
energies = collect(range(-3.1, length=500, stop=3.1))
Γ = 0.01
DOS = Spectrum.dos_dense(getbloch(hops), ks, energies; Γ=Γ)

plot(energies, DOS)
mkpath("output")
savefig("output/dos.pdf")