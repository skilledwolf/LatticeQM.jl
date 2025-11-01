# Concept Guide — Geometry & Lattices

This guide expands on Tutorial 1 with additional context about how LatticeQM
represents lattices, handles dimensionality, and exposes helper utilities.

## Core types
- `Structure.Lattice`: stores primitive vectors (`A`), fractional coordinates,
  and optional extra dimensions (e.g. layer labels, z offsets).
- `Structure.Geometries`: collection of convenience constructors for common
  lattices (`honeycomb`, `triangular`, `square`, multilayer stacks).
- `Structure.Lattices`: lower-level factory supporting programmatic lattice
  generation.

## Constructing lattices
```julia
lat = Lattice([1 0; 0 1], [[0, 0] [0.5, 0.5]])   # square lattice with two sites
lat3d = Lattice([[1,0,0] [0,1,0] [0,0,1]])       # cubic lattice
twisted = Geometries.honeycomb_twisted(13)       # moiré honeycomb (θ ≈ 1.05°)
```
- Use the optional `extra_dimensions` keyword to attach metadata:
  `Lattice(A, coords; extra_dimensions=["layer"])`.
- Prefer `Geometries` helper functions for well-tested configurations; inspect
  their docstrings for supported keywords (strain, offsets, stacking).

## Neighbours and symmetry
- `Structure.getneighbors(lat; cutoff, order)` returns neighbour shells with
  displacement vectors—essential for building hoppings.
- `Structure.brillouinzone(lat)` and `Structure.highsymmetrypoints(lat)` help
  generate k-paths aligned with lattice symmetries.
- `Structure.symmetries(lat)` (where available) exposes point-group data;
  leverage it to reduce computational workloads.

## Visualisation
```julia
using Plots
plot(lat; repeat=[0:1, 0:1], bonds=true, legends=false)
```
- For multilayer systems, use the `LayeredLayouts` recipes in `Plotting` to
  separate layers vertically.
- Export figures alongside scripts (`mkpath("output"); savefig("output/lat.svg")`)
  to keep provenance.

## Best practices
- Avoid mutating lattice internals directly; instead use
  `Structure.update_spacecoordinates!` or `Structure.update_A!`.
- Document any custom lattice constructors and add them to `Geometries` when
  stable so users can discover them via the API reference.
- Cross-reference the `extra/examples/` scripts to showcase non-standard
  lattices (e.g. TMDs, triangular Hubbard models).
