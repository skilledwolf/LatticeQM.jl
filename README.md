# LatticeQM.jl

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17500140.svg)](https://doi.org/10.5281/zenodo.17500140) [![Docs: master](https://img.shields.io/badge/docs-master-blue.svg)](https://tobiaswolf.net/LatticeQM.jl/master/)

LatticeQM is a Julia package for lattice structures and tight-binding models: build operators, compute bands and expectation values, and study mean-field, linear response, and Floquet physics. Licensed under AGPL-3.0-only.

## Installation
```julia
using Pkg
Pkg.add("LatticeQM")
```

## Documentation
- See `extra/docs` for the documentation source. Build with:
  `julia --project=extra/docs -e 'include("extra/docs/make.jl")'`, then open `extra/docs/build/index.html`.

## Citation
If you use LatticeQM in research or publications, please cite:

```
Tobias Wolf (2025). skilledwolf/LatticeQM.jl: v0.8.0 (v0.8.0). Zenodo. https://doi.org/10.5281/zenodo.17500140
```

## Acknowledgements
- I thank Dr. Oded Zilberberg and Dr. Gianni Blatter for their guidance and support as my thesis advisors; I developed the package in the course of my PhD research. 
- Tobias Kästli (@vigoleis) contributed the Floquet module as part of his Master's project and authored the accompanying tutorial.

## Author
- Dr. Tobias Wolf
