# LatticeQM.jl

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17500140.svg)](https://doi.org/10.5281/zenodo.17500140) [![Docs: master](https://img.shields.io/badge/docs-master-blue.svg)](https://tobiaswolf.net/LatticeQM.jl/master/)

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
| `:multithreaded`   | `julia -t N` / `JULIA_NUM_THREADS` | dense diagonalization on a multi-core machine        |
| `:distributed`     | `julia -p N` / `addprocs(N)`       | sparse / Arpack workloads, or many cheap k-points    |

A few caveats worth knowing about:

- **BLAS oversubscription.** Each Julia worker (and each Julia thread, via libblastrampoline) opens its own BLAS thread pool. Running `addprocs(N)` on an N-core machine therefore launches roughly N×N BLAS threads, which thrashes the cache. The fix is to pin BLAS to one thread per Julia worker:
  ```julia
  using Distributed
  addprocs(N)
  using LinearAlgebra
  @everywhere using LinearAlgebra
  @everywhere LinearAlgebra.BLAS.set_num_threads(1)
  ```
  The same applies when using `multimode=:multithreaded` with a Julia thread count higher than 1: set BLAS to 1 thread to leave the linear-algebra parallelism to Julia.
- **`:multithreaded` requires a dense format.** Arpack (used for sparse diagonalization) is not thread-safe. Pass `format=:dense` or use `multimode=:distributed` for sparse problems.
- **Workers must `using LatticeQM` before being used.** Either `julia -p N --project=.` and `@everywhere using LatticeQM`, or `addprocs(N; exeflags=\"--project=$(Base.active_project())\")` followed by `@everywhere using LatticeQM`.

## Citation
If you use LatticeQM in research or publications, please cite:

```
Tobias Wolf (2025). skilledwolf/LatticeQM.jl: v0.8.0 (v0.8.0). Zenodo. https://doi.org/10.5281/zenodo.17500140
```

## Acknowledgements
- I thank Dr. Oded Zilberberg and Dr. Gianni Blatter for their guidance and support as my thesis advisors; I developed the package in the course of my PhD research. 
- Tobias Kästli (@vigoleis) contributed the Floquet module as part of his Master's project and authored the accompanying tutorial.

## Author
- Dr. Tobias Wolf
