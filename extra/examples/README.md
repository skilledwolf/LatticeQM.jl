# LatticeQM examples

Working scripts that exercise the public API. Most run in seconds on a
laptop; the SCF and sweep ones take minutes to hours and benefit from
distributed mode (`julia -p N`).

For polished tutorials see `extra/tutorial/` instead — this folder is the
research-style scratch space, kept runnable.

## Index

### `graphene/` — single-layer graphene
| File | What | Runtime | Parallelism |
|---|---|---|---|
| [`dos.jl`](graphene/dos.jl) | Density of states via `Spectrum.getdos` | sec | serial |
| [`hofstadter.jl`](graphene/hofstadter.jl) | Hofstadter butterfly via Peierls substitution | sec | serial |
| [`opticalconductivity.jl`](graphene/opticalconductivity.jl) | `LinearResponse.opticalconductivity` | sec | serial |
| [`ribbons.jl`](graphene/ribbons.jl) | Armchair/zigzag ribbons via `reduceto1D` | sec | serial |
| [`ribbon_haldane.jl`](graphene/ribbon_haldane.jl) | Haldane edge states on a ribbon | sec | serial |
| [`valleybands.jl`](graphene/valleybands.jl) | Valley-resolved bands | sec | serial |
| [`hubbardmeanfield_example1.jl`](graphene/hubbardmeanfield_example1.jl) | SCF Hubbard with Zeeman field | min | threaded |
| [`hubbardmeanfield_example2.jl`](graphene/hubbardmeanfield_example2.jl) | SCF Hubbard with sublattice imbalance | min | threaded |
| [`hubbardmeanfield_supercell.jl`](graphene/hubbardmeanfield_supercell.jl) | SCF on a holey supercell | many min | distributed |
| [`hubbardsweep.jl`](graphene/hubbardsweep.jl) | Sweep gap-vs-U via `progress_map` | many min | threaded |
| [`non_unitary_sc.jl`](graphene/non_unitary_sc.jl) | (U,V) BdG sweep | hour+ | threaded |
| [`non_unitary_sc_plot.jl`](graphene/non_unitary_sc_plot.jl) | Plot script for the above | sec | serial |

### `graphene_AB/` — AB-stacked bilayer graphene
The `system.jl` helper builds the standard Hamiltonian; `sampling.jl`
defines a refined-around-K Brillouin-zone sampler.
| File | What | Runtime |
|---|---|---|
| [`dos.jl`](graphene_AB/dos.jl) | DOS on a custom non-uniform k-grid | min |
| [`fermisurface.jl`](graphene_AB/fermisurface.jl) | Fermi surface scan, primitive cell | min |
| [`meanfield_single.jl`](graphene_AB/meanfield_single.jl) | Single SCF run, save JLD2 | min |
| [`fermisurfaces/main.jl`](graphene_AB/fermisurfaces/main.jl) | Three FS variants (graphene TB / rhombohedral × 2 windows), pick via `LATTICEQM_PARAMS` | min |
| [`sweep/meanfield_sweep.jl`](graphene_AB/sweep/meanfield_sweep.jl) | μ-sweep + SCF, six parameter sets in [`params.jl`](graphene_AB/sweep/params.jl) | hours |
| [`sweep/post_fermisurface.jl`](graphene_AB/sweep/post_fermisurface.jl) | FS post-processing on each converged ρ | min |

### `graphene_ABC/` — ABC-stacked trilayer graphene
| File | What |
|---|---|
| [`sampling.jl`](graphene_ABC/sampling.jl) | Custom K-point sampler |
| `*.ipynb` | Exploration notebooks |

### `triangular/`
| File | What | Runtime |
|---|---|---|
| [`hubbardmeanfield_example1.jl`](triangular/hubbardmeanfield_example1.jl) | 120° non-collinear order on a triangular lattice | min |

### `twistedgraphene/` — twisted bilayer (single-shot)
| File | What | Runtime |
|---|---|---|
| [`valleybands.jl`](twistedgraphene/valleybands.jl) | Valley-resolved bands at twist angle n=5 | min |
| [`valleybands_displace.jl`](twistedgraphene/valleybands_displace.jl) | Same with vertical displacement modulation | min |

### `twistedgraphene_scf/` — TBG SCF (laptop / few-core)
[`main.jl`](twistedgraphene_scf/main.jl) runs Hubbard SCF on twisted bilayer
graphene. Parameter sets are listed at the top of the file; pick one with
`julia main.jl <name>` or `LATTICEQM_PARAMS=<name>`. Runs in ~5 min for
the default `:testing` set on 8 workers.

### `twistedgraphene_slurm/` — TBG SCF on a SLURM cluster
Same calculation as `twistedgraphene_scf/` but tuned for production runs
on TACC Lonestar6. See [`README.md`](twistedgraphene_slurm/README.md) for
the cluster setup steps; [`submission_script.sh`](twistedgraphene_slurm/submission_script.sh)
is a working sbatch template.

### `SSHxSSH/`
| File | What |
|---|---|
| [`wilsonloop.jl`](SSHxSSH/wilsonloop.jl) | Nested Wilson loops on a 2D SSH model (Benalcazar–Bernevig–Hughes) |

## Conventions

- **Output directories** named `output_*` are gitignored. Examples create
  them on demand. Feel free to `rm -rf` between runs.
- **Parameter sets**: scripts that explore multiple parameter choices use a
  `PARAMS = Dict(...)` at the top and pick the active set via
  `julia <script>.jl <name>` (positional arg) or `LATTICEQM_PARAMS=<name>`
  (env var, more SLURM-friendly). The default if neither is given is
  whatever the script declares as its first / smallest variant.
- **`multimode`**: examples accept `:serial` / `:multithreaded` /
  `:distributed`; most pick one based on `nworkers()`. Run with
  `julia -p N --project=../../.. <script>.jl` for distributed.
- **`BLAS.set_num_threads(1)`**: required for distributed/threaded modes
  to avoid BLAS oversubscription on multi-core boxes. Already set in the
  examples that need it.
- **Apple Silicon / macOS performance**: the default OpenBLAS backend is
  significantly slower than Apple's Accelerate framework for the dense
  Hermitian `heevr` calls that dominate the chempot solve and density
  matrix on big TBG cells. To switch:

  ```julia
  ]add AppleAccelerate          # one-off
  using AppleAccelerate         # MUST come before `using LatticeQM`
  using LatticeQM
  ```

  On the n=10 TBG cell (2648 orbitals), this is a **~3× speedup on
  serial dense `eigen!`** and ~1.8× on threaded — the single biggest
  performance lever on M-series.

  If you're on OpenBLAS and stuck with it, the threaded path benefits
  from `BLAS.set_num_threads(3)` on a 9-thread Julia (`-t 9`) — about
  30% faster than BLAS=1. Pass via
  `Parallel.configure_blas!(exec; threads=3)` if you call it directly,
  or set `BLAS.set_num_threads(3)` and toggle
  `LatticeQM.Parallel._BLAS_CONFIGURED[] = true` to suppress the
  later default-pin to 1. Apple Accelerate doesn't benefit from
  BLAS>1 — its `heevr` doesn't internally thread.

## Smoke tests

The fast examples are exercised end-to-end by `test/test_examples_smoke.jl`,
gated on `LATTICEQM_FULL_EXAMPLES=1` so they don't slow down `Pkg.test()`
by default. To run them:

```sh
LATTICEQM_FULL_EXAMPLES=1 julia --project=. -e 'using Pkg; Pkg.test()'
```
