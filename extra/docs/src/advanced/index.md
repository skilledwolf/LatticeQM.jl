# Advanced Workflows

Guidance for scaling calculations, running on clusters, and organising outputs.

## Performance playbook
- Profile hotspots with Julia's built-in tools (`@time`, `@allocated`, `Profile`).
- Prefer sparse vs. dense representations according to system size; use
  `TightBinding.hopdim` to estimate memory needs before allocating.
- For very large problems, consider chunking k‑grids and streaming results to
  disk (HDF5/JLD2) for decoupled plotting.

## Parallel & HPC execution
- Twisted bilayer workloads benefit from HPC resources. Use
  `extra/examples/twistedgraphene_slurm/submission_script.sh` as a starting
  template, and update resource requests (nodes, walltime) once validated.
- Record successful job configurations (cluster name, Julia/BLAS settings) to
  help others reproduce results.

## Data management
- Store heavy outputs (HDF5, JLD2) under clearly named folders (e.g.
  `output/<date>_<description>/`). Keep repository footprint manageable by
  tracking derived artefacts through git-lfs or external storage where needed.
- Include README files in data directories summarising provenance and parameter
  choices.

## Troubleshooting checklist
- **Long runtimes**: Verify sparse vs. dense representations, reduce k-point
  density, or truncate Floquet harmonics.
- **Memory spikes**: Inspect intermediate allocations (use `@time`, `@allocated`);
  split calculations into batches if necessary.
- **Numerical instabilities**: Tighten or loosen convergence tolerances, switch
  mixing strategies, or seed from prior converged states.

## Work-in-progress tracking
- Track gaps and ideas in your issue tracker; link scripts and figures so
  others can reproduce and review changes.
