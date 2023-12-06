using Distributed
nprocs() < 1 && addprocs(3);

using LatticeQM

lat = Lattices.superlattice(Geometries.honeycomb_AB(), [[1,1] [2,-1]])
Lattices.foldPC!(lat; shift=[1/3,0,0])#; shift=[1/3,1/3,0])

H = Operators.graphene_rhombohedral(lat; γ0=-1.0, γ1=0.3, γ2=-0.045, γ3=0.07, γ4=0.045)
Operators.addinterlayerbias!(H, lat, 0.3)
# H = Operators.graphene(lat; tz=0.3, ℓinter=0.1, ℓintra=0.05, cellrange=3, format=:dense)
# Operators.addinterlayerbias!(H, lat, 0.3)

posZ = Operators.positionalong(lat, [0,0,1.0]; rescale=true)
valley = Operators.valley(lat; spinhalf=false)

# # or for square grid in real space
B = Lattices.getB(lat)[1:2,1:2]
kpoints = 0.1*(Structure.regulargrid(; nk=300^2) .- [0.5,0.5])
kgrid = inv(B)*kpoints

# Calculate bands on fine grid
bands, obs = Spectrum.bandmatrix(H, kgrid, [valley, posZ]; multimode=:distributed)

# Energies at which to calculate the fermi surfaces
energies_fs = [-0.045, -0.047, -0.055, -0.065]
# energies_fs = [-0.095, -0.107, -0.125, -0.2]
broadening = 0.0002#0.0003

density = Spectrum.fermisurfacedensity_fromdata(bands, energies_fs; broadening=broadening)
density_valley = Spectrum.fermisurfacedensity_fromdata(bands, energies_fs, obs[:,:,1]; broadening=broadening)
density_layer = Spectrum.fermisurfacedensity_fromdata(bands, energies_fs, obs[:,:,2]; broadening=broadening)#; broadening=0.001)

using DelimitedFiles
writedlm("kpoints.out", kpoints)
writedlm("fs_energy.out", energies_fs)
writedlm("fs_density.out", density)
writedlm("fs_density_layer.out", density_layer)
writedlm("fs_density_valley.out", density_valley)
