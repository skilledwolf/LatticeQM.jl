# Tutorials Overview

LatticeQM provides a sequence of Jupyter notebooks in `extra/tutorial/` that
walk through progressively more advanced workflows. Each documentation page in
this section mirrors a notebook, highlighting the learning goals, main API
calls, and checkpoints to verify that results look sensible.

## Running the notebooks
1. Launch Julia in the project environment:
   ```julia
   julia --project=.
   using Pkg; Pkg.instantiate()
   ```
2. Start Jupyter (or VS Code) with the instantiated environment:
   ```julia
   using IJulia; notebook(dir="extra/tutorial")
   ```
3. Execute cells in order. When a tutorial saves artefacts (plots, HDF5 files,
   logs), they appear alongside the notebook or within its `output/` subfolder.

If you prefer a headless run to validate the notebooks, you can automate
notebook execution with a simple Julia script and review the generated
artefacts.

## Tutorial roadmap
| # | Topic | Notebook | Highlights | Notes |
| --- | --- | --- | --- | --- |
| 1 | Structure & Geometry | `Tutorial1_Structure.ipynb` | Build lattices from scratch, explore predefined geometries, render 2D layouts. | Establishes terminology used throughout the docs. |
| 2 | Band Structures | `Tutorial2_Bands.ipynb` | Construct graphene Hamiltonians, sample Brillouin-zone paths, plot band dispersions. | Reuses geometries from Tutorial 1. |
| 3 | Haldane Model & Topology | `Tutorial3_Haldane.ipynb` | Add topological mass terms, compute Chern numbers, visualise edge modes. | Introduces topology utilities. |
| 4 | Twisted Bilayer Graphene | `Tutorial4_Twisted.ipynb` | Build large moiré lattices, compare continuum and lattice approaches, manage memory footprint. | `Tutorial4_Twisted2/3` capture experimental variations—treat their outputs as provisional until revalidated. |
| 5 | Hofstadter Butterfly | `Tutorial5_Hofstadter.ipynb` | Thread magnetic flux, generate Hofstadter spectra, analyse fractal band structures. | Heavy k-point sampling; expect longer runtimes. |
| 6 | Mean-field Self-Consistency | `Tutorial6_Meanfield.ipynb` | Run Hartree–Fock loops, inspect convergence metrics, export density matrices. | Shares utilities with `extra/examples/graphene/hubbardmeanfield_*`. |
| 7 | Floquet Dynamics | `Tutorial7_Floquet.ipynb` | Construct time-periodic drives, compute quasienergy spectra, monitor resonances. | Builds on modules co-authored with Tobias Kästli. |

For background reading and additional implementation details, cross-reference
the **Concept Guides** and the repository examples under `extra/examples/`.
