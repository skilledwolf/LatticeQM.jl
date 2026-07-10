# LatticeQM.jl

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17500139.svg)](https://doi.org/10.5281/zenodo.17500139) [![Docs: master](https://img.shields.io/badge/docs-master-blue.svg)](https://tobiaswolf.net/LatticeQM.jl/master/)

LatticeQM is a Julia package for lattice structures and tight-binding models: build operators, compute bands and expectation values, and study mean-field, linear response, and Floquet physics. Licensed under AGPL-3.0-only.

## Installation
```julia
using Pkg
Pkg.add("LatticeQM")
```

## Documentation
- See `extra/docs` for the documentation source. Build with:
  `julia --project=extra/docs -e 'include("extra/docs/make.jl")'`, then open `extra/docs/build/index.html`.

## Parallelism

Several routines (`bandmatrix`/`getbands`, `getdos`, `getdensitymatrix!`, …) accept a `multimode` keyword:

| `multimode`        | Activated by               | When to use                                                  |
|--------------------|----------------------------|--------------------------------------------------------------|
| `:serial`          | always available            | small grids, debugging                                       |
| `:multithreaded`   | `julia -t N` / `JULIA_NUM_THREADS` | dense or sparse diagonalization on a multi-core machine |
| `:distributed`     | `julia -p N` / `addprocs(N)`       | many k-points spread across worker processes         |
| `:auto` (default for most routines) | always available | picks `:distributed` if workers exist, else threads, else serial |

A few caveats worth knowing about:

- **BLAS oversubscription is handled automatically.** Each Julia worker (and each Julia thread, via libblastrampoline) opens its own BLAS thread pool, so `addprocs(N)` on an N-core machine would otherwise launch roughly N×N BLAS threads and thrash the cache. Parallel routines call `Parallel.configure_blas!` internally on first use, pinning BLAS to one thread per worker/thread for `:multithreaded` and `:distributed` runs (`:serial` is left untouched). To override, call `LatticeQM.Parallel.configure_blas!(exec; threads=k)` yourself before the first parallel call — the automatic configuration is idempotent and won't undo your choice.
- **Sparse + threaded is supported.** Sparse diagonalization uses KrylovKit, which is pure Julia and thread-safe, so `format=:sparse` works with every `multimode`, including `:multithreaded`.
- **Workers must `using LatticeQM` before being used.** Either `julia -p N --project=.` and `@everywhere using LatticeQM`, or `addprocs(N; exeflags=\"--project=$(Base.active_project())\")` followed by `@everywhere using LatticeQM`.

## Citation
If you use LatticeQM in research or publications, please cite:

```
Tobias Wolf. skilledwolf/LatticeQM.jl. Zenodo. https://doi.org/10.5281/zenodo.17500139
```

This concept DOI always resolves to the latest release; see the [Zenodo record](https://doi.org/10.5281/zenodo.17500139) for version-specific DOIs.

## Acknowledgements
- I thank Dr. Oded Zilberberg and Dr. Gianni Blatter for their guidance and support as my thesis advisors; I developed the package in the course of my PhD research. 
- Tobias Kästli (@vigoleis) contributed the Floquet module as part of his Master's project and authored the accompanying tutorial.

## Author
- Dr. Tobias Wolf
