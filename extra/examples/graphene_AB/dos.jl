
include("system.jl")

lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=false, spinhalf=false)

include("sampling.jl")
kgrid, kweights = getcustomsampling(; N=5, k=2, l=500)
energies = LinRange(-0.08,-0.17, 500)
dos = Spectrum.getdos(H, energies, kgrid, kweights; format=:dense, Γ=0.002)

using DelimitedFiles
mkpath("output_dos")
writedlm("output_dos/energies.out", energies)
writedlm("output_dos/dos.out", dos)