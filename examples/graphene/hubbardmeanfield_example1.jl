using LinearAlgebra, Plots
using LatticeQM
using LatticeQM.Operators: graphene, gethaldane, valleyoperator, addsublatticeimbalance!, addzeeman!
using LatticeQM.Operators: magnetization, density
using LatticeQM.Meanfield

lat = Geometries2D.honeycomb()
sx, sy, sz, sublA, sublB = getoperator(lat, ["SX", "SY", "SZ", "sublatticeAspin", "sublatticeBspin"])
ks = kpath(lat; num_points=200)

hops = graphene(lat; mode=:spinhalf)
addzeeman!(hops, lat, r->sign(r[4]-0.5).*1.5.*[sin(0.0π),0,cos(0.0π)] )

# Set up interaction
v = gethubbard(lat; mode=:σx, a=0.5, U=5.0) # interaction potential
ρ_init = initialguess(v, :random; lat=lat) # initial guess
hf = hartreefock(hops, v)

ρ_sol, ϵ_GS, HMF, converged, error = solveselfconsistent( # run the calculation
    hf, ρ_init, 0.75; klin=30, iterations=800, tol=1e-7,# p_norm=Inf,
    T=0.01, β=0.25,  show_trace=true, clear_trace=true
)

# Get magnetization
mA, mB = real.(magnetization(ρ_sol, [sublA,sublB], lat))
δM = mA - mB; M = mA+mB
Mabs = norm(M); δMabs = norm(δM)
dens = density(ρ_sol)
@info("Groundstate energy", ϵ_GS)
@info("Magnetization", Mabs, δMabs, M, δM)
@info("Density", dens)

# Get the bands without mean-field terms
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, sz)
p1 = plot(bands; markersize=2, size=(600,200))

# Get the bands with mean-field terms
bands_mf = getbands(HMF.h, ks, sz)
bands_mf.bands .-= HMF.μ # shift chemical potential to zero
p2 = plot(bands_mf; markersize=2, size=(600,200))

# Show the band structure side-by-side
plot!(p1, title="B=1.5, non-interacting")
plot!(p2, title="B=1.5, Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

mkpath("output"); savefig("output/hubbardmeanfield_example1.pdf")
