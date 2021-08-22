# using Distributed
# addprocs(5)

using LinearAlgebra
using SparseArrays
using LatticeQM

######################################################################
# Prepare system and output
######################################################################

Us = -LinRange(0,3.5,20)
Vs = LinRange(0,1.0,20)
filling = 0.48
klin = 50

######################################################################
# Run distributed calculation
######################################################################
using ProgressBars
println(":: Load data...")


using DelimitedFiles
BANDGAP = readdlm("output/non_unitary_sc/bandgap.out")
SWAVEAUP = readdlm("output/non_unitary_sc/swaveAup.out")
SWAVEBUP = readdlm("output/non_unitary_sc/swaveBup.out")
SWAVEADOWN = readdlm("output/non_unitary_sc/swaveAdown.out")
SWAVEBDOWN = readdlm("output/non_unitary_sc/swaveBdown.out")
DENSITY0 = readdlm("output/non_unitary_sc/density.out")
DENSITYdiff = readdlm("output/non_unitary_sc/density_imbalance.out")

######################################################################
# Create and save plot
######################################################################
println(":: Creating and saving plot...")

using Plots

p1 = heatmap(-Us,Vs, log.(BANDGAP'); clims=(-3,1), xlabel="-U", ylabel="V", title="Bandgap", size=(300,300))
p2 = heatmap(-Us,Vs,log.(abs.(DENSITYdiff')); clims=(-6,1), xlabel="-U", ylabel="V", title="Density A-B", size=(300,300))
p3 = heatmap(-Us,Vs, log.(SWAVEAUP+SWAVEADOWN+SWAVEBUP+SWAVEBDOWN)'; clims=(-6,1), xlabel="-U", ylabel="V", title="S-WAVE (A+B)-(UP+DOWN)", size=(300,300))
p4 = heatmap(-Us,Vs, log.(abs.(SWAVEAUP+SWAVEADOWN-(SWAVEBUP+SWAVEBDOWN)))'; clims=(-6,1), xlabel="-U", ylabel="V", title="S-WAVE |A-B|-(UP+DOWN)", size=(300,300))
# p5 = heatmap(-Us,Vs, log.(SWAVEAUP+SWAVEBUP)'; clims=(-4,0),  xlabel="-U", ylabel="V", title="S-WAVE (A+B)-UP", size=(300,300))
# p6 = heatmap(-Us,Vs,(SWAVEBUP-SWAVEBUP)'; xlabel="-U", ylabel="V", title="S-WAVE (A-B)-UP", size=(300,300))
# p7 = heatmap(-Us,Vs, (SWAVEADOWN+SWAVEBDOWN)'; xlabel="-U", ylabel="V", title="S-WAVE (A+B)-DOWN", size=(300,300))
# p8 = heatmap(-Us,Vs,(SWAVEBDOWN-SWAVEBDOWN)'; xlabel="-U", ylabel="V", title="S-WAVE (A-B)-DOWN", size=(300,300))

plot(p1,p2,p3,p4; layout=(2,2), size=(600,500))

mkpath("output")
savefig("output/non_unitary_sc.pdf")