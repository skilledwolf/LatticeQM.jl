
# using AppleAccelerate
using LoopVectorization
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1) # set number of threads for BLAS

using Plots
using LatticeQM
import LatticeQM.Structure.Lattices

# println("Number of workers: ", nworkers())

println("Generate lattice...")
lat = Geometries.honeycomb_twisted(5)
print("Number of sites: ", Lattices.countorbitals(lat), "\n")
Lattices.foldPC!(lat)
sx, sy, sz, sublA, sublB = Operators.getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])

println("Generate lattice Hamiltonian...")
hops = Operators.graphene(lat; format=:sparse, mode=:spinhalf, tz=0.52)
Operators.addzeeman!(hops, lat, 1e-5)
# Operators.addzeeman!(hops, lat, r->sign(r[4]-0.5).*1.5.*[sin(0.0π),0,cos(0.0π)] )

filling = 0.5 + 6.0 / hopdim(hops)
println("Target filling: ", filling)

# Set up interaction
println("Set up meanfield interaction...")
v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=2.7) # interaction potential
ρ_init = Meanfield.initialguess(v, :random, :Z; lat=lat) # initial guess


println("Create band plot...")
# Get the bands without mean-field terms
ks_plot = kpath(lat; num_points=180)

Operators.setfilling!(hops, filling; nk=9, multimode=:distributed)
bands_plot = getbands(hops, ks_plot, sz; format=:sparse, num_bands=36, multimode=:serial)

mkpath("output_n5")
plot(bands_plot; colorbar=true)
savefig("output_n5/hubbardmeanfield_example1_noninteracting.pdf")
# exit()


println("Run selfconsistent solver...")
ρ_sol, ϵ_GS, HMF, converged, residue = Meanfield.solvehartreefock( # run the calculation
    hops, v, ρ_init, filling; klin=6, iterations=25, tol=1e-3,# p_norm=Inf,
    T=0.002, β=0.93, show_trace=true, clear_trace=true, verbose=true, multimode=:distributed
)


println("Get magnetization...")
mA, mB = real.(Operators.magnetization(ρ_sol, [sublA,sublB], lat))
δM = mA - mB; M = mA+mB
Mabs = norm(M); δMabs = norm(δM)
dens = Operators.density(ρ_sol)
@info("Groundstate energy", ϵ_GS)
@info("Magnetization", Mabs, δMabs, M, δM)
@info("Density", dens)

println("Create band plot...")
# Get the bands without mean-field terms
ks = kpath(lat; num_points=170)

Operators.setfilling!(hops, filling; nk=9, multimode=:distributed)
bands = getbands(hops, ks, sz; format=:sparse, num_bands=36, multimode=:distributed)
p1 = plot(bands; markersize=1, size=(600,200), colorbar=true)

# Get the bands with mean-field terms
hmf = HMF.h
Operators.setfilling!(hmf, filling; nk=9, multimode=:distributed)
bands_mf = getbands(hmf, ks, sz; format=:sparse, num_bands=36, multimode=:distributed)
# bands_mf.bands .-= HMF.μ # shift chemical potential to zero
p2 = plot(bands_mf; markersize=1, size=(600,200), colorbar=true)

# Show the band structure side-by-side
plot!(p1, title="noninteracting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

mkpath("output_n5"); savefig("output_n5/hubbardmeanfield_example1_distributed.pdf")

println("Exiting...")
sleep(5)