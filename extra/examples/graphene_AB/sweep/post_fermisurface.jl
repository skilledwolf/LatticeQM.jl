using Distributed
nprocs() < 2 && addprocs(7)

using LinearAlgebra
BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "..", "system.jl"))
include(joinpath(@__DIR__, "params.jl"))

using FileIO, DelimitedFiles

println("Building system...")
lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=true, spinhalf=true)

# Read sweep results from `meanfield_sweep.jl`.
isdir(OUTDIR) || error("$OUTDIR not found — run meanfield_sweep.jl first.")
μs       = readdlm(joinpath(OUTDIR, "chemicalpotentials.out"))
fillings = readdlm(joinpath(OUTDIR, "fillings.out"))

# Square 0.1×0.1 patch in (b1,b2) BZ coordinates centered at Γ.
B = Lattices.getB(lat)[1:2,1:2]
kpoints = 0.1 * (Structure.regulargrid(; nk=300^2) .- [0.5,0.5])
kgrid   = inv(B) * kpoints

println("Starting calculation...")
for (i_, filling) in enumerate(fillings)
    println("  doing #$i_ out of $(length(fillings))")

    HMF = FileIO.load(joinpath(OUTDIR, "meanfield$i_.jld2"), "HMF")

    bands, obs = Spectrum.bandmatrix(HMF.h, kgrid, ops; multimode=:distributed)

    bandgrid_dir = joinpath(OUTDIR, "bandgrid_$i_")
    mkpath(bandgrid_dir)
    writedlm(joinpath(bandgrid_dir, "kpoints.out"),    kpoints)
    writedlm(joinpath(bandgrid_dir, "bandmatrix.out"), bands)
    writedlm(joinpath(bandgrid_dir, "obsmatrix.out"),  obs)

    energies = [HMF.μ]
    density        = Spectrum.fermisurfacedensity_fromdata(bands, energies;             broadening=0.0005)
    density_layer  = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,1]; broadening=0.0005)
    density_valley = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,2]; broadening=0.0005)
    density_spin   = Spectrum.fermisurfacedensity_fromdata(bands, energies, obs[:,:,3]; broadening=0.0005)

    fs_dir = joinpath(OUTDIR, "fs_$i_")
    mkpath(fs_dir)
    writedlm(joinpath(fs_dir, "kpoints.out"),       kpoints)
    writedlm(joinpath(fs_dir, "density.out"),       density)
    writedlm(joinpath(fs_dir, "density_layer.out"), density_layer)
    writedlm(joinpath(fs_dir, "density_valley.out"),density_valley)
    writedlm(joinpath(fs_dir, "density_spin.out"),  density_spin)
end

println("Done.")
