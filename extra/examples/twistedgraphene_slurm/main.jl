
using Distributed
# using AppleAccelerate

@everywhere begin
  using LinearAlgebra#, MKL
  LinearAlgebra.BLAS.set_num_threads(1) # set number of threads for BLAS
  using LatticeQM
end

using Plots
import LatticeQM.Structure.Lattices

# n_angle, tz, outname, U_hubbard, guess = 11, 0.46, "output_n11_fm", 2.5, :ferro
# n_angle, tz, outname, U_hubbard, guess = 6, 0.52, "output_testing", 2.6, :ferro
# n_angle, tz, outname, U_hubbard, guess = 10, 0.52, "output_n10_fm3", 2.2, :ferro
# n_angle, tz, outname, U_hubbard, guess = 10, 0.52, "output_n10_fm4", 2.8, :ferro
# n_angle, tz, outname, U_hubbard, guess = 13, 0.39, "output_n13_fm5", 2.8, :ferro
n_angle, tz, outname, U_hubbard, guess = 13, 0.39, "output_n13_random", 2.4, :random

multimode = (nworkers() > 1) ? :distributed : :multithreaded

@info "Generate lattice..."
lat = Geometries.honeycomb_twisted(n_angle)
println("Number of sites: ", Lattices.countorbitals(lat))
sx, sy, sz, sublA, sublB = Operators.getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])

@info "Generate Hamiltonian..."
hops = Operators.graphene(lat; format=:sparse, mode=:spinhalf, tz=tz)
Operators.addzeeman!(hops, lat, 1e-9)
# Operators.addzeeman!(hops, lat, r->sign(r[4]-0.5).*1.5.*[sin(0.0π),0,cos(0.0π)] )

filling = 0.5 + 6.0 / hopdim(hops)
println("Target filling: ", filling)

# Set up interaction
@info "Set up mean-field interaction..."
v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=U_hubbard) # interaction potential

# ρ_init = Meanfield.initialguess(v, :random, :Z; lat=lat) # initial guess
ρ_init = Meanfield.initialguess(v, guess; lat=lat) # initial guess

@info "Band plot..."
# Get the bands without mean-field terms
ks_plot = kpath(lat; num_points=100)

# @info "... setting filling"
@everywhere (hops = $hops, ks_plot = $ks_plot, sz = $sz, v=$v)
Operators.setfilling!(hops, filling; nk=9^2, multimode=multimode)

# @info "... getting bands (sparse)"
bands_plot = getbands(hops, ks_plot, sz; format=:sparse, num_bands=36, multimode=multimode)

mkpath("$outname")
plot(bands_plot; colorbar=true)
savefig("$outname/bands_nonint.pdf")
# exit()


@info "Run mean field..."
ρ_sol, ϵ_GS, HMF, converged, residue = Meanfield.solvehartreefock( # run the calculation
    hops, v, ρ_init, filling; klin=9, iterations=16, tol=5e-3,# p_norm=Inf,
    T=0.002, β=0.93, show_trace=true, verbose=false, multimode=multimode
)
@everywhere GC.gc()
# exit()

@info "Magnetization..."
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

Operators.setfilling!(hops, filling; nk=9^2, multimode=multimode)
bands = getbands(hops, ks, sz; format=:sparse, num_bands=36, multimode=multimode)
p1 = plot(bands; markersize=1, size=(600,200), colorbar=true)

# Get the bands with mean-field terms
hmf = Meanfield.hMF(HMF)
Operators.addchemicalpotential!(hmf, -HMF.μ)

@everywhere (hmf=$hmf)
bands_mf = getbands(hmf, ks, sz; format=:sparse, num_bands=36, multimode=multimode)
p2 = plot(bands_mf; markersize=1, size=(600,200), colorbar=true)

# Show the band structure side-by-side
plot!(p1, title="noninteracting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

mkpath("$outname");
savefig("$outname/bands_meanfield.pdf");

println("Exiting...")
# sleep(5)
