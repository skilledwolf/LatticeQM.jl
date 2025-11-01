# Tutorial 8 — Superconductivity (BdG)

Model spin-singlet superconductivity via self-consistent Bogoliubov–de Gennes
(BdG) mean-field on a single-layer graphene lattice. This example mirrors the
parameters from an interacting honeycomb study (onsite attraction with mild
nearest‑neighbour interaction, slight doping), and plots quasiparticle bands
with particle/hole branches and the pairing gap.

## Learning goals
- Build a spinful graphene Hamiltonian and add a weak sublattice mass.
- Construct a BdG operator and an attractive onsite interaction for s‑wave pairing.
- Run a self-consistent BdG (Hartree–Fock) cycle at slight doping.
- Plot the BdG bands and inspect the gap opening.

## Setup
```@setup sc
using LatticeQM, Plots
lat = Geometries.honeycomb()
# Spinful graphene tight-binding
h0 = Operators.graphene(lat; mode=:spinhalf, format=:dense)

# Build BdG operator embedding
HBdG = Superconductivity.BdGOperator(h0)

# Short-ranged interaction: onsite attraction (U<0) + mild NN term (V>0)
# Values adapted from calculation_interacting2.jl
V = Operators.getshortrangedpotential(lat, -2.0, 0.25; spin=true)

# Build a BdG density-matrix initial guess (random, nonlocal seed)
ρ0_init = Meanfield.initialguess(V, :random, :nonlocal; lat=lat)
Δ0_init = Meanfield.initialguess(V, :random, :nonlocal; lat=lat)
ρ_init = Superconductivity.BdGOperator(ρ0_init, Δ0_init)

# Slight electron doping (fraction of occupied states per spin)
filling = 0.49

# k-sampling for self-consistency (use legacy-like density for closer parity)
klin = 70
```

## Self-consistent BdG
```@example sc
ρ_sc, ϵ_GS, HF, converged, resid = Meanfield.solvehartreefock(
    HBdG, V, ρ_init, filling; klin=klin, iterations=500, tol=1e-5,
    T=0.001, β=0.75, show_trace=false
)
println("Converged: ", converged, ", residual = ", resid)
```

```@example sc
# Extract the full BdG mean-field Hamiltonian and chemical potential
HBdG_mf = Meanfield.hMF(HF)
μ = HF.μ
occ = Spectrum.filling(HBdG_mf, μ; nk=klin)
println("μ = ", round(μ, digits=4), ", filling ≈ ", round(occ, digits=3))
size(HBdG_mf([0.0, 0.0]))
```

## Bands and pairing gap
```@example sc
figdir = joinpath(pwd(), "figures"); mkpath(figdir)
ks = kpath(lat; num_points=200) # Γ–K–M–Γ path
# Shift by chemical potential for plotting parity with legacy script
Operators.addchemicalpotential!(HBdG_mf, -μ)
# Colour by electron weight using an electron-sector projector function
eP = k -> Superconductivity.electron(HBdG_mf)(k)
bands_bdg = Spectrum.getbands(HBdG_mf, ks, eP)
# Plot symmetric window around zero to highlight particle/hole branches
p = plot(bands_bdg, 1; marker=:none, size=(520, 320),
         ylabel="ε/t", title="Graphene BdG bands (Δ ≠ 0)",
         colorbar=true, colorbar_title="electron weight")
savefig(p, joinpath(figdir, "graphene_bdg_bands.svg"))
nothing
```

![](figures/graphene_bdg_bands.svg)

Tip: Increase |U| or k‑point density for a larger/cleaner gap at the cost of
runtime. To verify s‑wave symmetry, inspect the local pairing observable via
`Operators.localobservables(ρ_sc, lat)`.

## Common pitfalls
- If self-consistency oscillates, reduce the mixing parameter (e.g. `β=0.5`) or
  initialise with a uniform seed instead of a random one.
- The sign convention for doping follows electron filling; values below 0.5
  indicate electron doping in the spinful model.
