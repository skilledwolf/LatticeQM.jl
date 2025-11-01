# Getting started

This page helps you install LatticeQM, pick an execution environment, and run a
five-minute graphene example that exercises the core APIs.

## Choose your environment
- **Docker quickstart** (recommended for a zero-friction demo).
- **VS Code + Remote Containers** (integrated development environment).
- **Native Julia installation** (for long-lived local setups and HPC nodes).

Regardless of setup, the tutorials and examples share the same project
environment, so data written inside containers is available on the host.

## Option A — Docker quickstart
1. Install [Docker](https://www.docker.com).
2. From the repository root run:
   ```bash
   docker-compose up --build
   ```
   This spins up a container with Julia, Python, and Jupyter configured.
3. Open the Jupyter URL printed in the terminal, navigate to `extra/tutorial`,
   and launch the desired notebook.

The folders `extra/examples` and `extra/tutorial` are mounted as volumes, so any
changes you make inside the container persist on the host filesystem.

## Option B — VS Code Remote Containers
1. Install Docker, [VS Code](https://code.visualstudio.com), and the
   *Remote - Containers* extension.
2. Open the repository in VS Code. The editor prompts you to reopen in a
   container; accept to reuse the same configuration as the Docker quickstart.
3. Use the integrated terminal or notebook UI to run tutorials and examples.

## Option C — Native Julia installation
1. Install the latest stable Julia release (≥1.10 recommended).
2. Add the package from the General registry:
   ```julia
   using Pkg
   Pkg.add("LatticeQM")
   ```
   For development work use `Pkg.develop(path="/path/to/LatticeQM")`.
3. Instantiate dependencies in the project environment:
   ```julia
   julia --project=.
   using Pkg; Pkg.instantiate()
   ```
4. Update periodically with `Pkg.update()`.

## Five-minute graphene example
Run the following snippet from a Julia REPL started with `--project=.` to verify
that everything is wired correctly.

```@example
using Plots
using LatticeQM

lat = Geometries.honeycomb()
hops = Operators.graphene(lat; mode=:spinhalf)
Operators.addzeeman!(hops, lat, 0.2)

valley = Operators.valley(lat; spinhalf=true)
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, valley)

plot(bands; ylabel="ε/t", colorbar_title="valley", size=(330, 200),
     colorbar=true, markercolor=:PiYG)
mkpath("figures") # hide
savefig("figures/graphene_bands_zeeman.svg") # hide
```

The script plots graphene bands with Zeeman splitting and stores the figure at
`figures/graphene_bands_zeeman.svg`.

## Troubleshooting
- **Package not found**: ensure you launched Julia with `--project=.` or called
  `Pkg.activate(".")`.
- **Plots backend errors**: install a compatible GR or Plotly backend
  (`Pkg.add("GR")`) or run inside the Docker container where it is preinstalled.
- **Long instantiation times**: precompile packages with
  `julia --project=. -e 'using Pkg; Pkg.precompile()'` before running notebooks.

## Next steps
- Continue with [Tutorial 1](tutorials/structure.md) for a deep dive into
  lattice construction.
- Explore scripts under `extra/examples/` for batch-ready workflows, or
  continue with additional tutorials.
