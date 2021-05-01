using Plots
using LatticeQM

# Load a lattice geometry
lat = Structure.Geometries.honeycomb()

# Construct graphene tight-binding Hamiltonian
hops = Operators.graphene(lat; mode=:spinhalf)

# Modify graphene Hamiltonian
#Operators.addhaldane!(hops, lat, 0.1; spinhalf=true)
Operators.addzeeman!(hops, lat, 0.2)

# Construct valley operator
valley = Operators.valleyoperator(lat; spinhalf=true)

# Get bands along (default) high-symmetry path
# and compute the expectation value of the valley operator for each eigenstate.
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, valley)

# Show bands
# save(bands, "bands.h5")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,200), colorbar=true, markercolor=:PiYG)

savefig("example.png")
