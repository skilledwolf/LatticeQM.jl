using LinearAlgebra, Plots
using LatticeQM

println("Generate lattice...")
lat = Structure.Geometries.honeycomb_twisted(6)
Structure.Lattices.foldPC!(lat)
sx, sy, sz, sublA, sublB = Operators.getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])
ks = kpath(lat; num_points=80)

println("Generate lattice Hamiltonian...")
hops = Operators.graphene(lat; format=:sparse, mode=:spinhalf, tz=0.45)
# Operators.addzeeman!(hops, lat, r->sign(r[4]-0.5).*1.5.*[sin(0.0π),0,cos(0.0π)] )

# Set up interaction
println("Set up meanfield interaction...")
v = Meanfield.gethubbard(lat; mode=:σx, a=0.5, U=2.0) # interaction potential
ρ_init = Meanfield.initialguess(v, :random; lat=lat) # initial guess

filling = 0.5+2.0/hopdim(hops)

println("Run selfconsistent solver...")
@time begin
    ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock( # run the calculation
        hops, v, ρ_init, filling; klin=2, iterations=100, tol=1e-6,# p_norm=Inf,
        T=0.002, β=0.93, show_trace=true, clear_trace=true, verbose=true, multimode=:distributed
    )
end

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
ks = kpath(lat; num_points=150)
Operators.setfilling!(hops, filling; nk=9, multimode=:distributed)
bands = getbands(hops, ks, sz; format=:sparse, num_bands=36, multimode=:distributed)
p1 = plot(bands; markersize=2, size=(600,200))

# Get the bands with mean-field terms
hmf = HMF.h
Operators.setfilling!(hmf, filling; nk=9, multimode=:distributed)
bands_mf = getbands(hmf, ks, sz; format=:sparse, num_bands=36, multimode=:distributed)
# bands_mf.bands .-= HMF.μ # shift chemical potential to zero
p2 = plot(bands_mf; markersize=2, size=(600,200))

# Show the band structure side-by-side
plot!(p1, title="noninteracting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

mkpath("output"); savefig("output/hubbardmeanfield_example1_distributed.pdf")
