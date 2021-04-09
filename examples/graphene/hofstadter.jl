using LatticeQM

#### Prepare system
lat = Geometries2D.honeycomb()
hops = Operators.graphene(lat; mode=:nospin, tz=0.3)
# Operators.addsublatticeimbalance!(hops, lat, 0.3)

#### Do the calculation
fluxes, energies = Operators.hofstadter(hops, lat, 40);

#### Produce plot
using PyPlot; const plt = PyPlot

points = hcat(([ϕ, ϵ] for (ϕ,en)=zip(fluxes,energies) for ϵ=en)...)

plt.scatter(points[1,:], points[2,:]; s=0.75)
plt.xlabel("flux per unit cell")
plt.ylabel("Energy")

mkpath("output")
plt.savefig("output/hofstadter.pdf")
plt.show()
println("Done!")