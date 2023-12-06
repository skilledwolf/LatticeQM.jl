
include("system.jl")

lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=true, spinhalf=false)

# kgrid = Structure.regulargrid(; nk=500^2)
# Lattices.foldBZ!(kgrid, lat)
# kpoints = Lattices.getB(lat) * kgrid
B = Lattices.getB(lat)[1:2,1:2]
kpoints = 0.1*(Structure.regulargrid(; nk=300^2) .- [0.5,0.5])
kgrid = inv(B)*kpoints

# Get band data
bands, obs = Spectrum.bandmatrix(H, kgrid, ops; multimode=:distributed2)


using DelimitedFiles
mkpath("output_fs2")
writedlm("output_fs2/kpoints.out", kpoints)
writedlm("output_fs2/bandmatrix.out", bands)
writedlm("output_fs2/obsmatrix.out", obs)

# Compute Fermi surfaces
energies = LinRange(-0.09,-0.16,20)
density = Spectrum.fermisurfacedensity_fromdata(bands, energies; broadening=0.0005)
density_layer = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,1]; broadening=0.0005)#; broadening=0.001)
density_valley = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,2]; broadening=0.0005)#; broadening=0.001)

writedlm("output_fs2/density.out", density)
writedlm("output_fs2/density_layer.out", density_layer)
writedlm("output_fs2/density_valley.out", density_valley)