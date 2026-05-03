# Bilayer-graphene Hubbard sweep

Runs a chemical-potential sweep over the AB-stacked bilayer graphene model
(`buildsystem` from `../system.jl`), solving the Hartree–Fock equations
self-consistently at each filling, then post-processes Fermi surfaces from
the converged mean-field Hamiltonians.

Replaces six near-duplicate folders that used to live under
`graphene_AB_sweep{1..5,1_reverse}/`. The variants are now selected via
`params.jl`; see that file for the full list.

## Run

```sh
# Mean-field sweep (writes output_<name>/meanfield<i>.jld2):
julia -p 7 --project=../../.. meanfield_sweep.jl sweep1

# Same, with a different parameter set:
julia -p 7 --project=../../.. meanfield_sweep.jl sweep4

# Or pick the variant via env:
LATTICEQM_PARAMS=sweep2 julia -p 7 --project=../../.. meanfield_sweep.jl

# Post-process Fermi surfaces (reads output_<name>/, writes the same dir):
julia -p 7 --project=../../.. post_fermisurface.jl sweep1
```

## Available parameter sets (in `params.jl`)

| Name | a | U | klin | iter | β | nk | warm-start | direction |
|---|---|---|---|---|---|---|---|---|
| `sweep1`  | 5.0 | 3.3 | 50  | 500  | 0.45 | 100 | yes | forward |
| `sweep1r` | 5.0 | 3.3 | 50  | 500  | 0.45 | 100 | yes | reverse |
| `sweep2`  | 4.0 | 3.6 | 50  | 900  | 0.45 | 100 | yes | forward |
| `sweep3`  | 5.0 | 3.3 | 75  | 900  | 0.35 | 100 | yes | forward |
| `sweep4`  | 5.0 | 3.3 | 120 | 1000 | 0.35 | 120 | no  | forward |
| `sweep5`  | 3.0 | 3.5 | 75  | 900  | 0.35 | 100 | yes | forward |

`a`, `U`: capped-Yukawa range and strength.
`klin`: linear k-grid for the SCF density-matrix integral.
`iter`/`β`: SCF iteration cap and Anderson mixing parameter.
`nk`: grid for `setfilling!`.
`warm-start`: reuse converged ρ from previous μ as initial guess for the next.
`direction`: walk fillings forward or reverse — different SCF basins near
phase transitions.
