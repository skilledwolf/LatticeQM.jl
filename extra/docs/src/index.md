# LatticeQM Documentation

LatticeQM is a Julia toolkit for building lattice models, constructing
tight-binding Hamiltonians, and analysing electronic structure across a broad
set of workflows (band structures, topology, mean-field, Floquet dynamics).

## Choose your path
- [Getting Started](getting_started.md) — installation, environment setup, and a
  five-minute graphene example.
- [Tutorials](tutorials/index.md) — guided Jupyter notebooks spanning lattice
  construction to Floquet physics.
  
- [Concept Guides](guides/geometry.md) — deeper discussions of geometry,
  operator construction, observables, and interactions.
- [Advanced Workflows](advanced/index.md) — performance tuning, batch runs,
  and data management strategies.
- [API Reference](api.md) — auto-generated listings of exported names once
  docstrings are standardised.
- [Contributor Guide](contributing.md) — development workflow, testing, and
  release coordination.

## What you'll find here
- **Physics-driven examples**: graphene, multilayer systems, twisted bilayers,
  Hubbard mean-field, and Floquet drives.
- **Executable walkthroughs**: tutorial pages run during the doc build,
  regenerating lightweight figures and outputs.
- **Extensible architecture**: modular subpackages (`Structure`, `Spectrum`,
  `Meanfield`, `Floquet`, etc.) designed to be recombined.

## Staying up to date
- When new capabilities land (e.g. superconductivity utilities, advanced Berry
  diagnostics), open an issue or discussion to suggest additional examples.

```@contents
Pages = [
    "getting_started.md",
    "tutorials/index.md",
    "guides/geometry.md",
    "advanced/index.md",
    "contributing.md",
    "api.md"
]
Depth = 2
```

## Acknowledgements
- I thank Dr. Oded Zilberberg and Dr. Gianni Blatter for their guidance and support as my thesis advisors; I developed the package in the course of my PhD research. 
- Tobias Kästli (@vigoleis) contributed the Floquet module as part of his Master's project and authored the accompanying tutorial.

## Author
- Dr. Tobias Wolf
