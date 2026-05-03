using Distributed
nprocs() < 2 && addprocs(3)

using LatticeQM

# Three parameter sets, replacing the old fermisurfaces, fermisurfaces_2,
# fermisurfaces_3 directories. The first uses the simple `graphene` tight-
# binding; the second/third use the rhombohedral (γ0..γ4) parametrisation
# with different bias and energy windows. Pick via
# `julia main.jl <name>` or `LATTICEQM_PARAMS=<name>`.
const PARAMS = Dict(
    :graphene_tb => (
        # Simple graphene_AB tight-binding; broad energy window.
        H_kind=:graphene,
        # graphene parameters:
        tz=0.3, ℓinter=0.1, ℓintra=0.05, cellrange=3,
        # mean-field bias:
        bias=0.3,
        energies_fs=[-0.095, -0.107, -0.125, -0.2],
        broadening=0.0002,
    ),
    :rhombohedral_high => (
        # Rhombohedral, large γ values, low-energy window.
        H_kind=:rhombohedral,
        γ0=-1.0, γ1=0.3, γ2=-0.045, γ3=0.07, γ4=0.045,
        bias=0.3,
        energies_fs=[-0.045, -0.047, -0.055, -0.065],
        broadening=0.0002,
    ),
    :rhombohedral_mid => (
        # Same rhombohedral parameters, mid-energy window.
        H_kind=:rhombohedral,
        γ0=-1.0, γ1=0.3, γ2=-0.045, γ3=0.07, γ4=0.045,
        bias=0.3,
        energies_fs=[-0.115, -0.13, -0.17],
        broadening=0.0002,
    ),
)

const params_name = Symbol(get(ENV, "LATTICEQM_PARAMS",
                                isempty(ARGS) ? "graphene_tb" : ARGS[1]))
haskey(PARAMS, params_name) ||
    error("Unknown param set: $params_name. Available: $(collect(keys(PARAMS)))")
const p = PARAMS[params_name]
const OUTDIR = "output_$(params_name)"
@info "Running param set" params_name p

lat = Lattices.superlattice(Geometries.honeycomb_AB(), [[1,1] [2,-1]])
Lattices.foldPC!(lat; shift=[1/3,0,0])

H = if p.H_kind === :graphene
    Operators.graphene(lat; tz=p.tz, ℓinter=p.ℓinter, ℓintra=p.ℓintra,
                            cellrange=p.cellrange, format=:dense)
elseif p.H_kind === :rhombohedral
    Operators.graphene_rhombohedral(lat; γ0=p.γ0, γ1=p.γ1, γ2=p.γ2,
                                          γ3=p.γ3, γ4=p.γ4)
else
    error("Unknown H_kind=$(p.H_kind)")
end
Operators.addinterlayerbias!(H, lat, p.bias)

posZ   = Operators.positionalong(lat, [0,0,1.0]; rescale=true)
valley = Operators.valley(lat; spinhalf=false)

# Square 0.1×0.1 patch in (b1,b2) BZ coordinates centered at Γ.
B = Lattices.getB(lat)[1:2,1:2]
kpoints = 0.1 * (Structure.regulargrid(; nk=300^2) .- [0.5,0.5])
kgrid   = inv(B) * kpoints

bands, obs = Spectrum.bandmatrix(H, kgrid, [valley, posZ]; multimode=:distributed)

density        = Spectrum.fermisurfacedensity_fromdata(bands, p.energies_fs;             broadening=p.broadening)
density_valley = Spectrum.fermisurfacedensity_fromdata(bands, p.energies_fs, obs[:,:,1]; broadening=p.broadening)
density_layer  = Spectrum.fermisurfacedensity_fromdata(bands, p.energies_fs, obs[:,:,2]; broadening=p.broadening)

using DelimitedFiles
mkpath(OUTDIR)
writedlm(joinpath(OUTDIR, "kpoints.out"),         kpoints)
writedlm(joinpath(OUTDIR, "fs_energy.out"),       p.energies_fs)
writedlm(joinpath(OUTDIR, "fs_density.out"),      density)
writedlm(joinpath(OUTDIR, "fs_density_layer.out"),density_layer)
writedlm(joinpath(OUTDIR, "fs_density_valley.out"),density_valley)

println("Wrote $OUTDIR/")
