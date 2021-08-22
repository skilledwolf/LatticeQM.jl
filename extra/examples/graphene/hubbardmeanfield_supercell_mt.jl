using LinearAlgebra, Plots
BLAS.set_num_threads(1)

import FileIO
using DelimitedFiles
using LatticeQM

function honeycombholes(; N=4,rad=0.2, rhombic=false)
    lat = Geometries.honeycomb()
    Lattices.translate!(lat, [1/3,1/3,0])
    lat = Lattices.superlattice(lat, [[N,0] [0,N]])
    Lattices.foldPC!(lat)#; shift=[1/3,-1/3,0])
    filterind = map(x->norm(x)>rad*norm(Lattices.getA(lat)[:,1]), eachcol(Lattices.positions(lat))) |> findall
    lat.spacecoordinates = lat.spacecoordinates[:,filterind]
    lat.extracoordinates = lat.extracoordinates[:,filterind]

    if rhombic
        Lattices.translate!(lat, [-1/2,-1/2,0])
        lat.spacecoordinates = mod.(lat.spacecoordinates, 1.0)
    end

    lat
end

@info "Building system"
lat = honeycombholes(N=5, rad=0.20, rhombic=true) # N=5, rad=0.2, # N=4, rad=0.2
sz = Operators.getoperator(lat,"sz")

println("Number of atoms: ", Lattices.countorbitals(lat))

@time hops = Operators.graphene(lat; mode=:spinhalf) #|> DenseHops
# Operators.addzeeman!(hops, lat, r->sign(r[4]-0.5).*1.5.*[sin(0.0π),0,cos(0.0π)] )
# Operators.addzeeman!(hops, lat, 1e-6)

@info "Mean field calculation"
# Set up interaction
v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=2.0) #|> DenseHops # interaction potential 
ρ_init = Meanfield.initialguess(v, :random; lat=lat) #|> DenseHops # initial guess 
# ρ_init = Hops()
# Operators.addzeeman!(ρ_init, lat, r->sign(r[4]-0.5)*[0,0,1] )

# ρ_init = FileIO.load("output/meanfield4.jld2", "ρ_sol")

@time ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock( # run the calculation
    hops, v, ρ_init, 0.5; klin=10, iterations=700, tol=1e-7,# p_norm=Inf,
    T=0.001, β=0.5,  show_trace=true, clear_trace=true, multimode=:multithread
)

mkpath("output")
FileIO.save("output/meanfield6.jld2", "ρ_sol", ρ_sol)

# # Get magnetization
M = real.(Operators.localmagnetization(ρ_sol, lat))
writedlm("output/M.out", M)
writedlm("output/XYZ.out", Lattices.positions(lat))

@info "Band structure plots"
ks = kpath(lat; num_points=200)
@time bands = getbands(hops, ks, sz)
@time bands_mf = getbands(HMF.h, ks, sz)
bands_mf.bands .-= HMF.μ # shift chemical potential to zero

Plots.default(show=false)
ENV["GKSwstype"]="nul"

p1 = plot(bands; markersize=2, size=(600,200))
p2 = plot(bands_mf; markersize=2, size=(600,200))

# Show the band structure side-by-side
plot!(p1, title="non-interacting")
plot!(p2, title="Hubbard meanfield")
plot(p1,p2, titlefont=font(8))

mkpath("output"); savefig("output/hubbardmeanfield_supercell.pdf")
