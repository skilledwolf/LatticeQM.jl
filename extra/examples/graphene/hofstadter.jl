using LatticeQM

#### Prepare system
lat = Geometries.honeycomb()
hops = Operators.graphene(lat; mode=:nospin, tz=0.3)
# Operators.addsublatticeimbalance!(hops, lat, 0.3)

#### Do the calculation
fluxes, energies = Operators.hofstadter(hops, lat, 40);

#### Produce plot
using Plots

points = hcat(([ϕ, ϵ] for (ϕ,en)=zip(fluxes,energies) for ϵ=en)...)

scatter(points[1,:], points[2,:]; markersize=0.75,
        xlabel="flux per unit cell", ylabel="Energy", legend=false)

mkpath("output")
savefig("output/hofstadter.pdf")
println("Done!")