# Concept Guide — Operators & Hamiltonians

The `LatticeQM.Operators` and `LatticeQM.TightBinding` modules provide tools to
build and manipulate tight-binding Hamiltonians. This guide complements
Tutorials 2–4 and the repository examples under `extra/examples/`.

## Core abstractions
- `TightBinding.Hops`: sparse dictionary mapping lattice displacements to
  hopping matrices. Use `DenseHops`/`SparseHops` wrappers to convert storage.
- `Operators.graphene(lat; mode=:spinhalf)`: convenience constructor for
  graphene-like models. Other helpers generate Rashba, Zeeman, or Haldane terms.
- `Operators.gethops(lat, f)`: generate hoppings from a lattice and a distance
  function `f(r₁, r₂)`.

## Building Hamiltonians
```julia
lat = Geometries.honeycomb()
hops = Operators.graphene(lat; mode=:spinhalf)
Operators.addzeeman!(hops, lat, 0.2)
Operators.addhaldane!(hops, lat, 0.1; spinhalf=true)
```
- Compose multiple terms by calling mutation helpers in sequence.
- Use `addhops!(Hops(), lat, f)` when building models from scratch; the helper
  automatically populates displacement keys.
- Convert to matrices with `TightBinding.getbloch(hops)` or pass `hops` to
  spectrum utilities directly.

## Observables & projectors
- Construct observables (valley, layer, spin) with dedicated helpers, e.g.
  `Operators.valley(lat; spinhalf=true)`.
- For density matrices, rely on utilities in `src/modules/Operators/densitymatrix.jl`.
  Document any advanced use (partial traces, symmetrisation) as you extend coverage.

## Managing sparsity
- Prefer sparse storage for large systems (`SparseHops(hops)`).
- Use `TightBinding.hopdim(hops)` to inspect matrix dimensions before allocating.
- Benchmark conversions for your problem sizes and record the settings that
  balance accuracy and runtime.

## Quality checks
- Verify Hermiticity by confirming `hops[δ]` and `hops[-δ]` are conjugate
  transposes.
- For custom terms, write smoke tests that diagonalise a representative system
  and compare to analytical expectations.

## Where to go next
- Combine with the Mean-field guide to incorporate interaction-driven terms.
- Reference the API page for the full list of exported constructors once docstrings
  are standardised.
