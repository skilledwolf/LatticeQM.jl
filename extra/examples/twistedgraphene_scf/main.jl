using Distributed
using AppleAccelerate
using LinearAlgebra
# LinearAlgebra.BLAS.set_num_threads(1)

using Plots
@everywhere using LatticeQM
import LatticeQM.Structure.Lattices

# Parameter sets: pick one via `julia main.jl <name>` or `LATTICEQM_PARAMS=<name>`.
# Default if no argument is `:testing`.
const PARAMS = Dict(
    :testing         => (n_angle=6,  tz=0.52, U=2.8, guess=:random, outname="output_testing"),
    :testing_fm         => (n_angle=6,  tz=0.52, U=2.8, guess=:ferro, outname="output_testing2"),
    :testing_purif   => (n_angle=6,  tz=0.52, U=2.8, guess=:ferro,  outname="output_testing_purif"),
    :n10_fm3         => (n_angle=10, tz=0.52, U=2.2, guess=:ferro,  outname="output_n10_fm3"),
    :n10_fm4         => (n_angle=10, tz=0.52, U=2.7, guess=:ferro,  outname="output_n10_fm4"),
    :n11_fm          => (n_angle=11, tz=0.46, U=2.5, guess=:ferro,  outname="output_n11_fm"),
)

const params_name = Symbol(get(ENV, "LATTICEQM_PARAMS",
                                isempty(ARGS) ? "testing" : ARGS[1]))
haskey(PARAMS, params_name) ||
    error("Unknown param set: $params_name. Available: $(collect(keys(PARAMS)))")
const p = PARAMS[params_name]
@info "Running param set" params_name p

multimode = (nworkers() > 1) ? :distributed : :multithreaded

@info "Generate lattice..."
lat = Geometries.honeycomb_twisted(p.n_angle)
println("Number of sites: ", Lattices.countorbitals(lat))
sx, sy, sz, sublA, sublB = Operators.getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])

@info "Generate Hamiltonian..."
hops = Operators.graphene(lat; format=:sparse, mode=:spinhalf, tz=p.tz)
Operators.addzeeman!(hops, lat, 1e-9)

filling = 0.5 - 6.0 / hopdim(hops)
println("Target filling: ", filling)

@info "Set up mean-field interaction..."
v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=p.U)
ρ_init = Meanfield.initialguess(v, p.guess; lat=lat) 

@info "Band plot..."
ks_plot = kpath(lat; num_points=100)

@everywhere (hops = $hops, ks_plot = $ks_plot, sz = $sz)
Operators.setfilling!(hops, filling; nk=9^2, multimode=multimode)
bands_plot = getbands(hops, ks_plot, sz; format=:sparse, num_bands=36, multimode=multimode)

mkpath(p.outname)
plot(bands_plot; colorbar=true)
savefig("$(p.outname)/bands_noninteracting.pdf")

@info "Run mean field..."
@everywhere (hops = $hops, v = $v)
ρ_sol, ϵ_GS, HMF, converged, residue = Meanfield.solvehartreefock(
    hops, v, ρ_init, filling; klin=6, iterations=30, tol=5e-3, hidebar=true,
    T=0.002, β=0.85, show_trace=true, verbose=false, multimode=multimode
)

@info "Magnetization..."
mA, mB = real.(Operators.magnetization(ρ_sol, [sublA,sublB], lat))
δM = mA - mB; M = mA+mB
Mabs = norm(M); δMabs = norm(δM)
dens = Operators.density(ρ_sol)
@info("Groundstate energy", ϵ_GS)
@info("Magnetization", Mabs, δMabs, M, δM)
@info("Density", dens)

println("Create band plot...")
ks = kpath(lat; num_points=170)

Operators.setfilling!(hops, filling; nk=9, multimode=multimode)
bands = getbands(hops, ks, sz; format=:sparse, num_bands=36, multimode=multimode)
p1 = plot(bands; markersize=1, size=(600,200), colorbar=true)

# Get the bands with mean-field terms
hmf = Meanfield.hMF(HMF)
Operators.addchemicalpotential!(hmf, -HMF.μ)

@everywhere (hmf=$hmf)
bands_mf = getbands(hmf, ks, sz; format=:sparse, num_bands=36, multimode=multimode)
p2 = plot(bands_mf; markersize=1, size=(600,200), colorbar=true)

plot!(p1, title="noninteracting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

savefig("$(p.outname)/bands_meanfield.pdf");

println("Exiting...")
