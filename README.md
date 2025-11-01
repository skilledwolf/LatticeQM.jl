# LatticeQM.jl

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

LatticeQM is a Julia package for lattice structures and tight-binding models: build operators, compute bands and expectation values, and study mean-field, linear response, and Floquet physics. Licensed under AGPL-3.0-only.

## Installation
```julia
using Pkg
Pkg.add(url="https://gitlab.com/skilledwolf/LatticeQM.jl.git")
```

## Documentation
- See `extra/docs` for the documentation source. Build with:
  `julia --project=extra/docs -e 'include("extra/docs/make.jl")'`, then open `extra/docs/build/index.html`.

## Acknowledgements
- I thank Dr. Oded Zilberberg and Dr. Gianni Blatter for their guidance and support as my thesis advisors; I developed the package in the course of my PhD research. 
- Tobias Kästli (@vigoleis) contributed the Floquet module as part of his Master's project and authored the accompanying tutorial.
IMPORTANT: This GitLab mirror is frozen (read‑only).

Active development has moved to GitHub:
https://github.com/skilledwolf/LatticeQM.jl

Please open issues and pull requests on GitHub. This repository remains as an
archived mirror for reference.

---
