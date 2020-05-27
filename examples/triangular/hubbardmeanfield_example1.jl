using LinearAlgebra, Plots
gr() # use the GR backend for plots
using LatticeQM

###################################################################################################
###################################################################################################

using LatticeQM.Operators: nearestneighbor!
using LatticeQM.TightBinding: addspin

# Set up lattice
lat = Geometries2D.triangular_supercell()
sz = getoperator(lat, "sz")
sub1, sub2, sub3 = [getoperator(lat, "sublattice", i, 2) for i=1:3]

# Get nearest-neighbor hops in the honeycomb lattice
hops = addspin(nearestneighbor!(Hops(), lat), :spinhalf)

###################################################################################################
###################################################################################################

using LatticeQM.Meanfield
using LatticeQM.Operators: magnetization

# Set up interaction
v = gethubbard(lat; mode=:σx, a=0.5, U=6.0) # interaction potential
ρ_init = initialguess(v, :random; lat=lat) # initial guess
hf = hartreefock(hops, v)

ρ_sol, ϵ_GS, HMF, converged, error = solveselfconsistent( # run the calculation
    hf, ρ_init, 0.5; klin=30, iterations=800, tol=1e-7,# p_norm=Inf,
    T=0.01, β=0.9,  show_trace=true, clear_trace=true
)

m1,m2,m3 = real.(magnetization(ρ_sol, [sub1,sub2,sub3], lat))
m = m1+m2+m3
Mabs = norm(m)
dens = Operators.density(ρ_sol)
@info("Groundstate energy", ϵ_GS)
@info("Magnetizations", Mabs, norm(m1), norm(m2), norm(m3))
@info("Magnetization vectors", Mabs, m, m1, m2, m3)
@info("Density", dens)

# Get the bands without mean-field terms
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, sz)
p1 = plot(bands; markersize=2, size=(600,200))

# Get the bands with mean-field terms
bands_mf = getbands(HMF.h, ks, sz)
bands_mf.bands .-= HMF.μ
p2 = plot(bands_mf; markersize=2, size=(600,200))

# Show the band structure side-by-side
plot!(p1, title="non-interacting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))
mkpath("scf_example1"); savefig("scf_example1/bandsmeanfield.pdf")

###################################################################################################
###################################################################################################

using DelimitedFiles
XYZ = transpose(Structure.positions(lat))
M = transpose(hcat([m1,m2,m3]...))
mkpath("scf_example1"); writedlm("scf_example1/positions.out", XYZ)
mkpath("scf_example1"); writedlm("scf_example1/magnetization.out", M)
