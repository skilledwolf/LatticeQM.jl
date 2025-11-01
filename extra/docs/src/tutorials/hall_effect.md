# Tutorial 3 — Haldane Model & Topology

**Notebook**: `extra/tutorial/Tutorial3_Haldane.ipynb`

Explore how to induce topological gaps in graphene-like systems and compute
observables that diagnose Chern insulating phases.

## Learning goals
- Introduce complex next-nearest-neighbour hoppings via `Operators.addhaldane!`.
- Study phase diagrams by sweeping Haldane mass and onsite potentials.
- Evaluate Berry curvature and Chern numbers with `Spectrum` utilities.
- Visualise chiral edge states on ribbon geometries.

## Prerequisites
- Outcomes from Tutorial 2 (band-structure workflows).
- Optional: knowledge of topological band theory terminology.

## Workflow outline
1. **Model setup** — Build a honeycomb lattice and baseline graphene
   Hamiltonian.
2. **Haldane term** — Apply `Operators.addhaldane!(hops, lat, φ;
   spinhalf=true)` and optionally stack with sublattice imbalance.
3. **Parameter scans** — Use simple loops or `range` objects to inspect phase
   transitions. Store observables with `@save` or by writing CSV/JSON.
4. **Topological invariants** — Invoke `Spectrum.chern(hops, ks; ...)` or
   related helpers to compute Chern numbers.
5. **Ribbon perspective** — Transform to ribbon geometries (see
   `extra/examples/graphene/ribbon_haldane.jl`) and plot edge-localised states.

## Live example
```@setup haldane
using LatticeQM, Plots, Statistics
lat = Geometries.honeycomb()
hops = Operators.graphene(lat; mode=:nospin, format=:dense)
Operators.addhaldane!(hops, lat, 0.05)
valley = Operators.valley(lat; spinhalf=false)
ks = kpath(lat; num_points=200)
```

```@example haldane
figdir = joinpath(pwd(), "figures")
mkpath(figdir)
nothing
```

```@example haldane
cherns = Spectrum.getcherns(hops, 15, 15)
println("Chern numbers (lowest two) = ", cherns[1:2])
```

```@example haldane
p = plot(
    getbands(hops, ks, [valley]);
    ylabel="ε/t",
    marker=:none,
    size=(420, 260),
    title="Haldane model band structure",
    colorbar=true,
    colorbar_title="valley"
)
savefig(p, joinpath(figdir, "haldane_bands.svg"))
nothing
```

![](figures/haldane_bands.svg)

## Berry curvature map
```@example haldane
kgrid, berry = Spectrum.berry(hops, 20, 20, [1])
vals = filter(!isnan, abs.(berry[:]))
maxval = isempty(vals) ? 1.0 : Statistics.quantile(vals, 0.98)
p = scatter(
    (Structure.Lattices.getB(lat) * hcat(reshape(kgrid, :, 1)...))[1, :],
    (Structure.Lattices.getB(lat) * hcat(reshape(kgrid, :, 1)...))[2, :];
    marker_z = vec(berry),
    markerstrokewidth=0,
    markersize=3,
    colorbar=true,
    clims=(-maxval, maxval),
    markercolor=:RdBu,
    xlabel="k₁",
    ylabel="k₂",
    size=(360, 280),
    legend=false,
    aspect_ratio=:equal
)
savefig(p, joinpath(figdir, "haldane_berry.svg"))
nothing
```

![](figures/haldane_berry.svg)

## Phase diagram (coarse grid)
```@example haldane
using LaTeXStrings
phi_vals = LinRange(-π, π, 25)
m_vals = LinRange(-1.2, 1.2, 25)
h0 = Operators.graphene(lat; mode=:nospin, format=:dense)
C = fill(0.0, length(phi_vals), length(m_vals))
for (i, ϕ) in enumerate(phi_vals), (j, Δ) in enumerate(m_vals)
    local_model = deepcopy(h0)
    Operators.addsublatticeimbalance!(local_model, lat, Δ)
    Operators.addhaldane!(local_model, lat, 0.2; ϕ=ϕ)
    C[i, j] = sum(Spectrum.berry(local_model, 10, 10, [1])[2])
end
p = heatmap(
    phi_vals ./ π,
    m_vals,
    transpose(round.(C; digits=2));
    xlabel=L"\phi/\pi",
    ylabel=L"m/t_2",
    colorbar_title="Chern",
    size=(420, 240),
    color=:seismic
)
savefig(p, joinpath(figdir, "haldane_phase.svg"))
nothing
```

![](figures/haldane_phase.svg)

## Ribbon dispersion
```@example haldane
N = 20
lat_ribbon = Structure.Lattices.reduceto1D(Geometries.honeycomb(), [[1, 1] [N, -N]])
h_ribbon = Operators.graphene(lat_ribbon; mode=:nospin, format=:dense, cellrange=1)
Operators.addhaldane!(h_ribbon, lat_ribbon, 0.25)
position = Operators.positionalong(lat_ribbon, Structure.Lattices.basis(lat_ribbon, 2); rescale=true, center=true)
ks_ribbon = kpath(lat_ribbon; num_points=140)
bands_ribbon = getbands(h_ribbon, ks_ribbon, position)
p = plot(
    bands_ribbon, 1;
    ylabel="ε/t",
    colorbar_title="transverse position",
    csymmetric=true,
    markersize=1.5,
    size=(450, 260)
)
savefig(p, joinpath(figdir, "haldane_ribbon.svg"))
nothing
```

![](figures/haldane_ribbon.svg)

## Nested Wilson Loops (SSH×SSH snapshot)
```@example haldane
using LatticeQM.Structure.Lattices: Lattice, addbasis!, addorbital!, addextra!
function SSHxSSH(t1X, t2X, t1Y, t2Y)
    lat = Lattice()
    addbasis!(lat, [1, 0]); addbasis!(lat, [0, 1]); addextra!(lat, "sublattice")
    addorbital!(lat, [0,   0,   2]); addorbital!(lat, [1/2, 0,   4])
    addorbital!(lat, [0, 1/2,   3]); addorbital!(lat, [1/2, 1/2, 1])
    lat.specialpoints = LatticeQM.Geometries.kdict_sq
    h = DenseHops()
    h[[0,0]] = zeros(ComplexF64, 4, 4)
    h[[0,0]][3,1] = t1X; h[[0,0]][2,3] = -t1Y; h[[0,0]][4,2] = t1X; h[[0,0]][1,4] = t1Y
    h[[0,0]] += h[[0,0]]'
    h[[1,0]] = zeros(ComplexF64, 4, 4); h[[1,0]][3,1] = t2X; h[[1,0]][2,4] = t2X; h[[-1,0]] = h[[1,0]]'
    h[[0,1]] = zeros(ComplexF64, 4, 4); h[[0,1]][2,3] = -t2Y; h[[0,1]][4,1] = t2Y; h[[0,-1]] = h[[0,1]]'
    lat, h
end
latSSH, hSSH = SSHxSSH(0.5, 1.0, 0.5, 1.0)
pol1, U1, pol2, U2 = Spectrum.NestedWilson2D(hSSH, 60, 60, 1:2)
p = plot(scatter(collect(1:length(pol1[:,1])), pol1[:,1]; ms=3.0, label="1"),
         scatter(collect(1:length(pol1[:,2])), pol1[:,2]; ms=3.0, label="2");
         size=(420, 260), xlabel="k₂ index", ylabel="polarisation k₁")
savefig(p, joinpath(figdir, "sshxssh_wilson.svg"))
nothing
```

![](figures/sshxssh_wilson.svg)

## Validation checklist
- Confirm Chern numbers match the expected ±1 quantisation when phases wrap.
- Create plots of Berry curvature hotspots and verify symmetry.
- Compare ribbon spectra with bulk band gaps for consistency.

## Suggested extensions
- Couple the model to linear-response routines in `LinearResponse` to compute
  Hall conductivities.
- Introduce disorder or electric fields to study robustness of edge channels.
- Export topological markers for later visualisation using scripts in
  `extra/examples/`.
