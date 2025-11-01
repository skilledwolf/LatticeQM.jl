# Concept Guide — Floquet Workflows

Periodic driving introduces additional structure beyond static lattice
Hamiltonians. The `Floquet` module (developed with Tobias Kästli during his Master's project) provides
convenience APIs for constructing and analysing driven systems.

## Core functions
- `Floquet.makefloquet(hops; harmonics, ω, gauge)`: build the truncated Floquet
  Hamiltonian by stacking harmonic sectors.
- `Floquet.getspectrum(floquet_hops; kwargs...)`: obtain quasienergies and
  states.
- `Floquet.makeobservable(op; harmonics)`: lift static observables to Floquet
  space for expectation-value calculations.

## Workflow checkpoints
1. **Static baseline** — Assemble the non-driven Hamiltonian using `Operators`.
2. **Drive specification** — Encode vector potentials or onsite modulations with
   helper functions (e.g. `Floquet.circular_drive`).
3. **Harmonic truncation** — Choose harmonic count based on drive strength;
   verify convergence by increasing the cutoff until spectra stabilise.
4. **Diagonalisation** — Solve for quasienergies; classify them modulo the drive
   frequency.
5. **Observables** — Evaluate micromotion, time-averaged expectation values, or
   Floquet-broadened densities of states.

## Numerical tips
- Use sparse representations; Floquet Hamiltonians are block-structured and can
  grow rapidly with harmonic count.
- Normalise phases consistently across harmonics to avoid discontinuities in
  observables.
- Profile memory/time with representative parameters; record findings in the
  Advanced Workflows section.

## Further reading
- Tutorial: [Floquet Dynamics](../tutorials/floquet.md) for a guided example.
- Examples: Extend the graphene scripts by adding periodic drives and compare
  static vs. driven observables.
  
