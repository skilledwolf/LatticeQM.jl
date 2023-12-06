
include("system.jl")

lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=false, spinhalf=false)

kgrid = Structure.regulargrid(; nk=300^2)
# Lattices.foldBZ!(kgrid, lat)
kpoints = Lattices.getB(lat) * kgrid

# Get band data
bands, obs = Spectrum.bandmatrix(H, kgrid, ops; multimode=:distributed2)


using DelimitedFiles
mkpath("output_fs")
writedlm("output_fs/kpoints.out", kpoints)
writedlm("output_fs/bandmatrix.out", bands)
writedlm("output_fs/obsmatrix.out", obs)

# Compute Fermi surfaces
energies = LinRange(-0.09,-0.16,20)
density = Spectrum.fermisurfacedensity_fromdata(bands, energies; broadening=0.004)
density_layer = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,1]; broadening=0.001)#; broadening=0.001)

writedlm("output_fs/density.out", density)
writedlm("output_fs/density_layer.out", density_layer)