# using Pkg; Pkg.activate(ENV["JuliaDev"]*"/LatticeQM/")
using LinearAlgebra, Plots
gr() # use the GR backend for plots
# using Revise
using LatticeQM

# Set up lattice
lat = Geometries2D.honeycomb()

# Get nearest-neighbor hops in the honeycomb lattice
using LatticeQM.Operators: graphene, addsublatticeimbalance!
hops = graphene(lat; mode=:nospin) # or mode=:spinhalf for spin-1/2
addsublatticeimbalance!(hops, lat, 0.1)

# Get band structure
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks)

plot(bands, size=(330,240)) # plot data
savefig("bands_oc.pdf")

using LatticeQM.LinearResponse: opticalconductivity
frequencies = collect(range(0.0, length=800, 5.0))
oc = opticalconductivity(frequencies, 1, 1, hops, lat; klin=500, T=0.015, Γ=0.02)
oc = round.(oc; digits=8)

# Plot the optical conductivity
p = plot(frequencies, real(oc) * π , label="Re")
plot!(p, frequencies, imag(oc) * π , label="Im")
mkpath("output"); savefig("output/opticalconductivity.pdf")
