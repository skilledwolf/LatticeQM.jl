using LatticeQM

#### Prepare system
lat = Geometries2D.honeycomb()
hops = Operators.graphene(lat; mode=:nospin, tz=0.3)
# Operators.addsublatticeimbalance!(hops, lat, 0.3)

#### Do the calculation
fluxes, energies = Operators.hofstadter(hops, lat, 40);


#### Produce plot
using Plots

p = plot()
for (ϕ,Es)=zip(fluxes,energies)
    scatter!(p, repeat([ϕ],length(Es)), Es; markersize=0.75, markercolor=:black, legend=false)
end
plot!(p, xlabel="flux per unit cell", ylabel="Energy E/t")
plot!(p, size=(400,300))

mkpath("output"); savefig(p, "output/hofstadter.pdf")
println("Done!")