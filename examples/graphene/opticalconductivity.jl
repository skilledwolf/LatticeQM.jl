using LinearAlgebra, Plots

# using Revise
using LatticeQM

# Set up lattice
lat = Geometries.honeycomb()

# Get nearest-neighbor hops in the honeycomb lattice
hops = Operators.graphene(lat; format=:dense, mode=:nospin) # or mode=:spinhalf for spin-1/2
hops = DenseHops(hops)
Operators.addsublatticeimbalance!(hops, lat, 0.001)

# Get band structure
ks = kpath(lat; num_points=200)
bands = Spectrum.getbands(hops, ks)

plot(bands, size=(330,240)) # plot data
mkpath("output"); savefig("output/bands_oc.pdf")

frequencies = LinRange(0.0, 5.0, 500)
Γ=0.01
oc = LinearResponse.opticalconductivity(frequencies, 1, 1, hops, lat; klin=1000, T=0.001, Γ=Γ)
oc = -(oc.-oc[begin]) ./ (frequencies .+ 1im*Γ) # Optical conductivity tensor

# Plot the optical conductivity
p = plot(frequencies, real(oc) * π , label="Re")
plot!(p, frequencies, imag(oc) * π , label="Im")
savefig("output/opticalconductivity.pdf")
