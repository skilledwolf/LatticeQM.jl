using LinearAlgebra, Plots
using LatticeQM
using LatticeQM.Operators: graphene, gethaldane, valleyoperator, addsublatticeimbalance!, addzeeman!
using LatticeQM.Operators: magnetization, density
using LatticeQM.Meanfield

println("Generate lattice...")
lat = Geometries2D.honeycomb_twisted(3)
sx, sy, sz, sublA, sublB = getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])
ks = kpath(lat; num_points=80)

println("Generate lattice Hamiltonian...")
hops = graphene(lat; format=:sparse, mode=:spinhalf, tz=0.45)
# addzeeman!(hops, lat, r->sign(r[4]-0.5).*1.5.*[sin(0.0π),0,cos(0.0π)] )

# Set up interaction
println("Set up meanfield interaction...")
v = gethubbard(lat; mode=:σx, a=0.5, U=5.0) # interaction potential
ρ_init = initialguess(v, :random; lat=lat) # initial guess
hf = hartreefock(hops, v)

println("Run selfconsistent solver...")
ρ_sol, ϵ_GS, HMF, converged, error = solveselfconsistent( # run the calculation
    hf, ρ_init, 0.5+2.0/hopdim(hops); klin=2, iterations=5, tol=1e-6,# p_norm=Inf,
    T=0.002, β=0.5, show_trace=true, clear_trace=false, verbose=true, parallel=true
)

println("Get magnetization...")
mA, mB = real.(magnetization(ρ_sol, [sublA,sublB], lat))
δM = mA - mB; M = mA+mB
Mabs = norm(M); δMabs = norm(δM)
dens = density(ρ_sol)
@info("Groundstate energy", ϵ_GS)
@info("Magnetization", Mabs, δMabs, M, δM)
@info("Density", dens)

println("Create band plot...")
# Get the bands without mean-field terms
ks = kpath(lat; num_points=80)
bands = getbands(hops, ks, sz; format=:sparse, num_bands=36)
p1 = plot(bands; markersize=2, size=(600,200))

# Get the bands with mean-field terms
bands_mf = getbands(HMF.h, ks, sz; format=:sparse, num_bands=36)
bands_mf.bands .-= HMF.μ # shift chemical potential to zero
p2 = plot(bands_mf; markersize=2, size=(600,200))

# Show the band structure side-by-side
plot!(p1, title="noninteracting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

mkpath("output"); savefig("output/hubbardmeanfield_example1.pdf")
