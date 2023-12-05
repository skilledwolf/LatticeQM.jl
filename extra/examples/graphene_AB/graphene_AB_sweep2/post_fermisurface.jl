using Distributed
nprocs() < 2 && addprocs(47)

include("system.jl")

using FileIO, DelimitedFiles

println("Building system...")
lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=true, spinhalf=true)

v = Operators.getcappedyukawa(lat; format=:dense, cellrange=3, spin=true, a=5.0, U=3.3) # interaction potential
ρ_init = Meanfield.initialguess(v, :random, :nonlocal; lat=lat) # initial guess

# Chemical potentials to sweep over
μs = readdlm("output_mf_sweep/chemicalpotentials.out")
fillings = readdlm("output_mf_sweep/fillings.out")

# Prepare k points for sampling
B = Lattices.getB(lat)[1:2,1:2]
kpoints = 0.1*(Structure.regulargrid(; nk=300^2) .- [0.5,0.5])
kgrid = inv(B)*kpoints

# Main calculations
println("Starting calculation...")
for (i_, filling) in enumerate(fillings)

    println("  doing #$i_ out of "*string(length(fillings)))

    HMF = FileIO.load("output_mf_sweep/meanfield$i_.jld2", "HMF")

    # Get band data
    bands, obs = Spectrum.bandmatrix(HMF.h, kgrid, ops; multimode=:distributed)

    using DelimitedFiles
    mkpath("output_bandgrid/$i_/")
    writedlm("output_bandgrid/$i_/kpoints.out", kpoints)
    writedlm("output_bandgrid/$i_/bandmatrix.out", bands)
    writedlm("output_bandgrid/$i_/obsmatrix.out", obs)

    # Compute Fermi surfaces
    energies = [HMF.μ]
    density = Spectrum.fermisurfacedensity_fromdata(bands, energies; broadening=0.0005)
    density_layer = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,1]; broadening=0.0005)#; broadening=0.001)
    density_valley = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,2]; broadening=0.0005)#; broadening=0.001)

    writedlm("output_fs/$i_/kpoints.out", kpoints)
    writedlm("output_fs/$i_/density.out", density)
    writedlm("output_fs/$i_/density_layer.out", density_layer)
    writedlm("output_fs/$i_/density_valley.out", density_valley)
end

println("Done.")