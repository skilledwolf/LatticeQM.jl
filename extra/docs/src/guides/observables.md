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

Select the occupied bands explicitly. For example, the lower SSH band has
centre 1/2 in the topological phase (intercell hopping 1, intracell hopping 0.4):

```julia
using LatticeQM, LinearAlgebra
H(k) = ComplexF64[0 0.4+cis(-2π*k[1]); 0.4+cis(2π*k[1]) 0]
centres, vectors, states = Spectrum.Wilson1D(H, 201, [1])
O = Diagonal([1.0, -1.0])
weighted = Spectrum.WilsonWeights1D(H, 201, O, [1])
# weighted.centres, weighted.weights, weighted.multiplicities
```

`wilsonweights` returns one entry per distinct Wilson eigenspace. At a
degenerate phase the weight is the total expectation of `O` over that
subspace, and `multiplicities` records its dimension. It does not assign
basis-dependent weights to individual degenerate bands. The Hamiltonian
wrapper counts repeated grid endpoints with half weight. When supplying
frames directly, pass `kgrid` for endpoint-inclusive or nonuniform grids:

```julia
kgrid, states = Spectrum.statesgrid1D(H, 201, [1])
weighted = Spectrum.wilsonweights(states, O; kgrid)
```

For spin/layer sectors separated by an operator gap at zero, use
`Spectrum.SectorWilson1D(H, 201, O, occupied_bands)`. It returns
`((negative=centres_minus, positive=centres_plus), Ogap)`. Same-sign
degeneracies are allowed. The projected operator gap must stay open; this
alone does not quantize the sector phases. Check the electronic gap and
momentum-grid convergence separately. The SSH sublattice operator used above
has a closed projected gap and therefore cannot define these sign sectors.

Nested polarizations use whole Wannier sectors. `WilsonWannier2D` chooses
phase cuts between the distinct phases at the first transverse momentum;
these must remain in gaps throughout the sampled loop. To choose wider
sectors, supply cuts for each nesting direction, for example:

```julia
p1, p2 = Spectrum.WilsonWannier2D(H2D, 41, 41, occupied_bands;
    wannier_cuts=([-0.25, 0.25], [-0.25, 0.25]))
```

These cuts separate phases near 0 from phases near 1/2 across the ±π branch
cut. Each output entry is the total nested polarization of one sector.
Crossings within a sector are allowed; crossings of its boundaries are
rejected when resolved by the grid. A winding nested polarization has no
single periodic average and is also rejected. Array inputs must omit the
repeated endpoints; the Hamiltonian wrapper removes them automatically.
The sector construction follows the
[nested Wilson-loop definition of Benalcazar, Bernevig, and Hughes](https://arxiv.org/abs/1611.07987).

All these routines assume orthonormal states in a fixed periodic orbital
basis. Compare phases modulo integers and unwrap continuous branches when
plotting. Quantization to 0 or 1/2 can be lost when the protecting
symmetry is broken. Inspect `extra/examples/SSHxSSH/wilsonloop.jl` for a
two-dimensional model.

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
