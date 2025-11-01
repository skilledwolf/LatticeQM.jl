# Concept Guide — Mean-field & Superconductivity

The `Meanfield` and `Superconductivity` modules implement Hartree–Fock and
Bogoliubov–de Gennes (BdG) workflows. This guide collects best practices and
links relevant examples.

## Mean-field toolkit
- `Meanfield.solveselfconsistent(problem; kwargs...)`: generic interface for
  iterative solvers (density, spin, pairing).
- `Meanfield.solvehartreefock(model; kwargs...)`: Hartree–Fock convenience
  wrapper returning a `HartreeFock` object.
- `Meanfield.initialguess(lat; mode=:random)` produces seed density matrices; a
  deterministic `:uniform` mode is available for benchmarking.
- `Meanfield.HartreeFock` struct exposes fields such as `density`, `energies`,
  `iterations`, and `converged`.

## Superconducting extension
- `Superconductivity.BdGOperator` builds Bogoliubov–de Gennes Hamiltonians
  from normal-state hoppings and pairing matrices.
- Pairing helpers (s-wave, d-wave) are available under
  `Superconductivity`—audit docstrings as you standardise them.

## Typical workflow
1. **Model definition** — Specify the lattice, base hopping matrices, and
   interaction parameters (onsite U, nearest-neighbour V, pairing strength).
2. **Initial guess** — Use `initialguess` or reuse converged densities from
   previous runs (`hf.density`).
3. **Iteration control** — Tune `maxiter`, `tol`, and mixing (linear,
   Anderson). Monitor `hf.residuals`.
4. **Post-processing** — Export densities, order parameters, and energies using
   `Operators.densitymatrix` utilities or custom writers.
5. **Visualisation** — Plot layer/sublattice resolved densities; cross-check
   with scripts under `extra/examples/` for consistency.

## Debugging convergence
- Inspect the residual history; oscillations often indicate the need for
  stronger damping or a better initial guess.
- Use `tmp_debug.jl` or dedicated sandbox notebooks to isolate problematic
  parameter sets before committing changes.
- Log solver metadata (mixing coefficients, runtime) to aid reproducibility.

## Resources
- Examples: `extra/examples/graphene/hubbardmeanfield_*`,
  `extra/examples/twistedgraphene_scf`.
- Tutorials: [Tutorial 6](../tutorials/meanfield.md) for a hands-on walkthrough.
