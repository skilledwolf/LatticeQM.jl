using LinearAlgebra, Plots
gr()

using LatticeQM

lat = Structure.Geometries.honeycomb()
H = Operators.graphene(lat; mode=:spinhalf)

ks = LatticeQM.Utils.regulargrid(nk=400^2)
energies = collect(range(-3.1, length=500, stop=3.1))
Γ = 0.01
DOS = Spectrum.getdos(H, ks, energies; format=:dense, Γ=Γ)

plot(energies, DOS)
mkpath("output")
savefig("output/dos.pdf")