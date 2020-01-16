# Imports
using LinearAlgebra, Plots
pyplot() # use the matplotlib backend for plots
using LatticeQM

# Set up lattice
lat = Geometries2D.honeycomb()

# Get nearest-neighbor hops in the honeycomb lattice
hops = Operatorsgraphene(lat; mode=:nospin) # or mode=:spinhalf for spin-1/2

# Compile Bloch Hamiltonian
h = getbloch(hops)

# Path in k-space
ks = Structure.kpath(lat; num_points=200)

# Get bandstructure
bands = getbands(h, ks)

save(bands, "bands_example.h5") # save raw data
p = plot(bands, size=(330,240)) # plot data
# display(p) # might need to be uncommented in non-interactive mode

savefig("bands_example.pdf") # save the figure
