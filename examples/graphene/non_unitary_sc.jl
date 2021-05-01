# using Distributed
# addprocs(5)

using LinearAlgebra
using SparseArrays
using LatticeQM

######################################################################
# Prepare system and output
######################################################################

println(":: Preparing system...")

lat = Geometries.honeycomb()
ks = kpath(lat; num_points=200)

hops = Operators.graphene(lat; mode=:spinhalf)

Us = -LinRange(0,3.5,20)
Vs = LinRange(0,1.0,20)
filling = 0.48
klin = 50

BANDGAP = zeros(length(Us),length(Vs))
SWAVEAUP = zeros(length(Us),length(Vs))
SWAVEBUP = zeros(length(Us),length(Vs))
SWAVEADOWN = zeros(length(Us),length(Vs))
SWAVEBDOWN = zeros(length(Us),length(Vs))
WAVEdiff = zeros(length(Us),length(Vs))
DENSITY0 = zeros(length(Us),length(Vs))
DENSITYdiff = zeros(length(Us),length(Vs))


######################################################################
# Run distributed calculation
######################################################################
using ProgressBars
println(":: Entering sweep...")

IJ = [(i_,j_) for i_=1:length(Us) for j_=1:length(Vs)]

lk = Threads.ReentrantLock()
Threads.@threads for (i_,j_)=ProgressBar(IJ)
    
    U = Us[i_]; V = Vs[j_]

    hopsBDG = BdGOperator(hops)
    v = Operators.getshortrangedpotential(lat, U, V)

    ρ0_init = initialguess(v, :random; kind=:nonlocal)
    Δ0_init = initialguess(v, :random; kind=:nonlocal)
    ρ_init = BdGOperator(ρ0_init, Δ0_init)

    # println("Starting solvehartreefock in $i_, $j_")

    ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock( # run the calculation
        hopsBDG, v, ρ_init, filling; klin=klin, iterations=1000, tol=1e-5,# p_norm=Inf,
        T=0.01, β=0.7,  show_trace=false, clear_trace=false
    )
    Operators.addchemicalpotential!(HMF.h, -HMF.μ)
    
    M, SC = Operators.localobservables(ρ_sol, lat)
    M = real(M)#round.(; digits=12)
    # SC = SC#round.(SC; digits=12)
    
    # save the data of interest
    # println("Writing result in $i_, $j_")
    lock(lk) do 
        BANDGAP[i_,j_] = Spectrum.bandgap(DenseHops(HMF.h.h); klin=30)
        SWAVEAUP[i_,j_] = abs(SC[1,1]+SC[4,1])
        SWAVEBUP[i_,j_] = abs(SC[1,2]+SC[4,2])
        SWAVEADOWN[i_,j_] = abs(SC[1,1]-SC[4,1])
        SWAVEBDOWN[i_,j_] = abs(SC[1,2]-SC[4,2])
        DENSITY0[i_,j_] = (M[1,1]+M[1,2])/2
        DENSITYdiff[i_,j_] = (M[1,1]-M[1,2])/2
    end

end


using DelimitedFiles
mkpath("output/non_unitary_sc/")
writedlm("output/non_unitary_sc/bandgap.out", BANDGAP)
writedlm("output/non_unitary_sc/swaveAup.out", SWAVEAUP)
writedlm("output/non_unitary_sc/swaveBup.out", SWAVEBUP)
writedlm("output/non_unitary_sc/swaveAdown.out", SWAVEADOWN)
writedlm("output/non_unitary_sc/swaveBdown.out", SWAVEBDOWN)
writedlm("output/non_unitary_sc/density.out", DENSITY0)
writedlm("output/non_unitary_sc/density_imbalance.out", DENSITYdiff)

######################################################################
# Create and save plot
######################################################################
println(":: Creating and saving plot...")

using Plots

p1 = heatmap(-Us,Vs, log.(BANDGAP'); clims=(-4,0), xlabel="-U", ylabel="V", title="Bandgap", size=(300,300))
p2 = heatmap(-Us,Vs,log.(abs.(DENSITYdiff')); clims=(-4,0), xlabel="-U", ylabel="V", title="Density A-B", size=(300,300))
p3 = heatmap(-Us,Vs, log.(SWAVEAUP+SWAVEADOWN+SWAVEBUP+SWAVEBDOWN)'; clims=(-4,0), xlabel="-U", ylabel="V", title="S-WAVE (A+B)-(UP+DOWN)", size=(300,300))
p4 = heatmap(-Us,Vs, log.(abs.(SWAVEAUP+SWAVEADOWN-(SWAVEBUP+SWAVEBDOWN)))'; clims=(-4,0), xlabel="-U", ylabel="V", title="S-WAVE |A-B|-(UP+DOWN)", size=(300,300))
p5 = heatmap(-Us,Vs, log.(SWAVEAUP+SWAVEBUP)'; clims=(-4,0),  xlabel="-U", ylabel="V", title="S-WAVE (A+B)-UP", size=(300,300))
p6 = heatmap(-Us,Vs,(SWAVEBUP-SWAVEBUP)'; xlabel="-U", ylabel="V", title="S-WAVE (A-B)-UP", size=(300,300))
p7 = heatmap(-Us,Vs, (SWAVEADOWN+SWAVEBDOWN)'; xlabel="-U", ylabel="V", title="S-WAVE (A+B)-DOWN", size=(300,300))
p8 = heatmap(-Us,Vs,(SWAVEBDOWN-SWAVEBDOWN)'; xlabel="-U", ylabel="V", title="S-WAVE (A-B)-DOWN", size=(300,300))

plot(p1,p2,p3,p4,p5,p6,p7,p8; layout=(4,2), size=(600,900))

mkpath("output")
savefig("output/non_unitary_sc.pdf")