# Tutorial 6 — Mean-field Self-Consistency

**Notebook**: `extra/tutorial/Tutorial6_Meanfield.ipynb`

Walk through self-consistent Hartree–Fock calculations on lattice models,
capturing density and order-parameter evolution.

## Learning goals
- Configure `Meanfield` problems using `solveselfconsistent` and
  `solvehartreefock`.
- Set up initial guesses with `Meanfield.initialguess` and customise convergence
  criteria.
- Monitor observables (charge imbalance, magnetisation) during iterations.
- Persist convergence logs and final density matrices for downstream analysis.

## Prerequisites
- Tutorials 1–2 for lattice and Hamiltonian construction.
- Optional: experience with the examples under
  `extra/examples/graphene/hubbardmeanfield_*`.

## Workflow outline
1. **Model definition** — Prepare the lattice and base Hamiltonian; define the
   interaction channels (e.g. Hubbard U, nearest-neighbour coupling).
2. **Initial state** — Use `initialguess(lat; mode=:random, seed=...)` or reuse
   previous converged densities.
3. **Solver invocation** — Call `solveselfconsistent(problem; maxiter, tol,
   mix=...)` capturing the returned `HartreeFock` object.
4. **Diagnostics** — Inspect `hf.iterations`, `hf.free_energy`, and custom logs
   written to disk; plot convergence trajectories.
5. **Post-processing** — Export observables and convert density matrices using
   `Operators.densitymatrix.jl` utilities if needed.

## Live example
```@setup meanfield
using LatticeQM, Plots, Random, LinearAlgebra, LatticeQM.Utils
Random.seed!(1)
lat = Geometries.honeycomb()
hops = Operators.graphene(lat; mode=:spinhalf)
# Observables (spin and sublattice projectors)
sx, sy, sz, sublA, sublB = Operators.getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])

# Non-interacting reference bands (coloured by sz)
ks = kpath(lat; num_points=200)
bands0 = getbands(hops, ks, sz)
p1 = plot(bands0; markersize=2, size=(380, 240), title="non-interacting")

# Set up interaction and initial guess
v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=4.0)
ρ_init = Meanfield.initialguess(v, :random; lat=lat)

# Solve mean-field (parameters aligned with the notebook, trimmed iterations)
ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock(
    hops, v, ρ_init, 0.499999;
    klin=60, iterations=300, tol=1e-5,
    T=0.01, β=0.7, show_trace=false
)
```

```@example meanfield
println("Converged: ", converged, " with residual ", error)
```

```@example meanfield
figdir = joinpath(pwd(), "figures")
mkpath(figdir)
nothing
```

```@example meanfield
sublA, sublB = Operators.getoperator(lat, ["sublatticeAspin", "sublatticeBspin"])
mA, mB = real.(Operators.magnetization(ρ_sol, [sublA, sublB], lat))
(mA, mB)
```

```@example meanfield
δM = mA - mB
M = mA + mB
println("|M| = ", norm(M), " |δM| = ", norm(δM))
```

```@example meanfield
# Mean-field bands (use HMF.hMF, shifted by μ), coloured by sz for comparison
bands_mf = getbands(HMF.hMF, ks, sz)
bands_mf.bands .-= HMF.μ
p2 = plot(bands_mf; markersize=2, size=(380, 240), title="Hubbard mean-field")
p = plot(p1, p2; size=(780, 260), titlefont=font(8))
savefig(p, joinpath(figdir, "meanfield_bands_compare.svg"))
nothing
```

![](figures/meanfield_bands_compare.svg)

## Triangular lattice snapshot
```@example meanfield
lat_tri = Geometries.triangular_supercell()
plot(lat_tri, size=(320, 220))
```

```@example meanfield
base_tri = Operators.nearestneighbor!(Hops(), lat_tri)
h_tri = Utils.dense(TightBinding.addspin(base_tri, :spinhalf))
v_tri = Operators.gethubbard(lat_tri; mode=:σx, a=0.5, U=4.0)
ρ_init_tri = Meanfield.initialguess(v_tri, :random; lat=lat_tri)
ρ_tri, _, HMF_tri, converged_tri, _ = Meanfield.solvehartreefock(
    h_tri, v_tri, ρ_init_tri, 0.5;
    klin=16, iterations=120, tol=1e-4, show_trace=false
)
println("Triangular converged: ", converged_tri)
```

```@example meanfield
ks_tri = kpath(lat_tri; num_points=140)
bands_tri = getbands(HMF_tri.hMF, ks_tri)
bands_tri.bands .-= HMF_tri.μ
p = plot(bands_tri; size=(360, 220), xlabel="k", ylabel="ε/t")
savefig(p, joinpath(figdir, "triangular_meanfield.svg"))
nothing
```

![](figures/triangular_meanfield.svg)

## Validation checklist
- Ensure convergence thresholds are satisfied and iterations do not oscillate.
- Compare final energies against reference runs (`extra/examples/graphene`).
- Save density matrices to `output/` and verify they can be reloaded.

## Suggested extensions
- Benchmark serial vs. threaded performance and note runtime/memory settings.
- Integrate with twisted bilayer calculations (`twistedgraphene_scf` example).
- Investigate superconducting order parameters by coupling to the
  `Superconductivity` module.
