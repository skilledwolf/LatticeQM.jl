# Tutorial 7 — Floquet Dynamics

**Notebook**: `extra/tutorial/Tutorial7_Floquet.ipynb`

Learn how to model periodically driven systems using the `LatticeQM.Floquet`
module and interpret quasienergy spectra.

## Learning goals
- Build Floquet Hamiltonians from static tight-binding models and driving
  protocols.
- Truncate harmonic spaces judiciously to capture resonant processes.
- Compute quasienergy spectra and track avoided crossings.
- Evaluate time-averaged observables and micromotion.

## Prerequisites
- Tutorials 1–2 for standard Hamiltonian workflows.
- Basic understanding of Floquet theory.

## Workflow outline
1. **Base model** — Start with a lattice and tight-binding Hamiltonian (graphene
   or a custom model).
2. **Drive specification** — Define periodic fields via helper constructors in
   `Floquet` (e.g. circularly polarised light, pulsed fields).
3. **Floquet construction** — Use `Floquet.makefloquet(hops; harmonics, ω)`
   or equivalent functions demonstrated in the notebook.
4. **Diagonalisation** — Solve for quasienergies with `Floquet.getspectrum` or
   the standard `Spectrum.getbands` applied to the extended Hamiltonian.
5. **Observables** — Analyse mode occupations, compute time-averaged currents,
   and export plots.

## Live example
```@setup floquet
using LatticeQM, Plots, LinearAlgebra
lat = Geometries.honeycomb()
hops = Operators.graphene(lat; mode=:spinhalf)
Hk(k) = Matrix(hops(k))
dim_h = size(Hk([0.0, 0.0]), 1)
drive = Floquet.periodicDrive(15.0)
Floquet.addcos!(drive, 0.08 .* Matrix{Float64}(I, dim_h, dim_h), 1)
HF = Floquet.FloquetOperator(k -> Hk(k), drive, 1)
ks = kpath(lat; num_points=36)
```

```@example floquet
Floquet.periodicDrive(5.0)
```

```@example floquet
figdir = joinpath(pwd(), "figures")
mkpath(figdir)
nothing
```

```@example floquet
bands = getbands(HF, ks)
Floquet.keepfirstFBZ!(bands, HF)
first(bands.bands, 3)
```

```@example floquet
p = plot(bands; ylabel="ε/t", marker=:none, size=(380, 240), title="Floquet quasienergies (first zone)")
savefig(p, joinpath(figdir, "floquet_bands.svg"))
nothing
```

![](figures/floquet_bands.svg)

## Chern numbers for driven graphene (toy amplitude)
```@example floquet
N = size(HF(ks.points[:, 1]), 1) ÷ 2
chern_vals = Spectrum.getcherns(HF, 16, 16, N-1:N+2)
chern_vals
```

## Two-band toy model
```@setup floquet_model
using LatticeQM, LinearAlgebra, LatticeQM.Utils, Plots
const sigma_pauli = LatticeQM.Utils.σs
μ = 1.0
a = 4.0
b = 1.5
J = 1.5
function HB(k; μ=μ, a=a, b=b, J=J)
    d1 = a * sin(2π * k[1])
    d2 = a * sin(2π * k[2])
    d3 = (μ - J) - 2b * (2 - cos(2π * k[1]) - cos(2π * k[2])) + J * cos(2π * k[1]) * cos(2π * k[2])
    sum(d .* sigma_pauli[i] for (i, d) in enumerate((d1, d2, d3)))
end
sqlat = Geometries.square()
path2 = kpath(sqlat; num_points=120)
bands_static = Spectrum.getbands(HB, path2)
```

```@example floquet_model
static_cherns = Spectrum.getcherns(HB, 12, 12)
static_cherns
```

```@example floquet_model
p = plot(bands_static; size=(360, 240), title="Static two-band model")
savefig(p, joinpath(pwd(), "figures", "floquet_static.svg"))
nothing
```

![](figures/floquet_static.svg)

## Driven two-band model
```@example floquet_model
ω = 12.0
D = 0.8
M = 3
drive2 = Floquet.periodicDrive(ω, [[D 0; 0 -D], [D 0; 0 -D]], [1, -1])
HF2 = Floquet.FloquetOperator(HB, drive2, M)
bands_driven = Spectrum.getbands(HF2, path2)
Floquet.keepfirstFBZ!(bands_driven, HF2)
p = plot(bands_driven; size=(360, 240), title="Driven bands", ylims=(-ω/2, ω/2))
savefig(p, joinpath(pwd(), "figures", "floquet_driven.svg"))
nothing
```

![](figures/floquet_driven.svg)

```@example floquet_model
chern_drive = Spectrum.getcherns(HF2, 12, 12, [2M + 1])
wind_drive = Spectrum.getwindnum(HF2, 12, 12, 2M + 1)
(chern_drive, wind_drive)
```

## Validation checklist
- Confirm quasienergy zones repeat modulo the drive frequency.
- Compare against static limits (`ω → ∞`) to ensure convergence.
- Save plots (typically under `output/floquet/`) and verify they match the
  notebook reference.

## Common pitfalls
- Insufficient harmonic truncation can create spurious gaps. Increase `M` until
  bands and observables stabilise.
- Beware energy folding: always compare spectra modulo `ω` and consider
  calling `Floquet.keepfirstFBZ!` before plotting.
- For large systems, prefer sparse types and limit the number of requested
  bands to keep memory in check.

## Suggested extensions
- Couple Floquet results to linear-response calculations to estimate pump-probe
  signatures.
- Explore different truncation strategies and document runtime vs. accuracy.
- Benchmark parallel execution where available and capture typical runtimes
  for your hardware and parameter choices.
