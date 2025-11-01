# Concept Guide — Observables & Spectral Analysis

Analyse bands, densities of states, Berry curvature, and linear-response
coefficients using the `Spectrum`, `LinearResponse`, and `Plotting` modules.

## Spectrum essentials
- `Spectrum.getbands(hops, ks, op=nothing)` computes eigenvalues and, if an
  operator is supplied, expectation values along the k-path.
- `Spectrum.getdos(hops, emin, emax; klin, Γ, format)` returns density-of-states
  estimates. For large systems use `format=:sparse`.
- `Spectrum.chern` and related routines derive topological invariants from
  eigenvectors or Wilson loops.

## Berry phases & Wilson loops
```julia
lat = Geometries.honeycomb()
hops = Operators.graphene(lat)
loops = Spectrum.wilsonloop(hops; nk=200, bands=1:2)
```
- Inspect the SSH×SSH example (`extra/examples/SSHxSSH/wilsonloop.jl`) for a
  working script.
- Normalise phases to maintain continuity; unwrap when plotting.

## Linear response
- `LinearResponse.opticalconductivity` and friends evaluate conductivities using
  Kubo formulas. See `extra/examples/graphene/opticalconductivity.jl`.
- Pay attention to broadening parameters and sum-rule checks.

## Plotting recipes
- `Plotting.plot(bands; kwargs...)` produces band plots with colour bars,
  annotated symmetry points, and legend customisation.
- For density-of-states, combine `Plots` with returned arrays to overlay
  multiple datasets.

## Best practices
- Cache k-paths and operators when scanning parameters to avoid recomputation.
- For extensive sweeps, stream results to disk (HDF5, JLD2) and perform plotting
  in a separate step to decouple compute from rendering.
- Document numerical tolerances (Γ, nk, integration step) alongside figures so
  others can reproduce them.
