# LatticeQM for Julia

This is a development repository for julia package that deals with lattice structures and tight-binding models.

Not intended for public use yet. Much of the documentation is still missing.  
To get started, check out the notebooks in the folder `tutorial`.  
If you want to know more, get in touch with me.


## Installation
```bash
$ julia
julia> ]
pkg> add <pkg_url>
```

or
```bash
$ julia
julia> using Pkg; Pkg.add("<pkg_url>")
```

You can either use the https protocol, in which case `<pkg_url>` should be `https://gitlab.ethz.ch/wolft/LatticeQM.jl.git`, or you can use
the SSH protocal, in which case `<pkg_url>` should be `git@gitlab.ethz.ch:wolft/LatticeQM.jl.git`. The first method is fool-proof, the second requires a properly set-up SSH key.

## Update

Do not forget to regularly pull updates of the package:
```bash
pkg> update
```

## Uninstalling

```bash
 pkg> remove <pkg_name>
```

## Contributing
You can either clone the repo yourself or use Julia's built-in convenience functions:
1. Run
   ```bash
   pkg> develop <pkg_url or name>
   ```
  
   (The files will then be located at `$JULIA_PKG_DEVDIR`, by default `~/.julia/dev/`)

2. Make your modifications, commit and push them

3. If you want to return to the repo version you can type
   ```bash
   pkg> free <pkg_url or name>
   ```


Tip: If you want to edit a specific module or function you can use, which opens the correct file:
```bash
julia> edit(myfunction) 
julia> edit(mymodule)
```

(The editor can be specified with the global variable `JULIA_EDITOR`)

(I haven't actually tested this myself. Let me know if there are problems.)

Some more info can be found [here](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) and [here](https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/11-developing-julia-packages).

## Example code

```julia
using Plots
using LatticeQM

# Load a lattice geometry
lat = Geometries2D.honeycomb()

# Construct graphene tight-binding Hamiltonian
hops = Operators.graphene(lat; mode=:spinhalf)

# Modify graphene Hamiltonian
#Operators.addhaldane!(hops, lat, 0.1; spinhalf=true)
Operators.addzeeman!(hops, lat, 0.2)

# Construct valley operator
valley = Operators.valleyoperator(lat; spinhalf=true)

# Get bands along (default) high-symmetry path
# and compute the expectation value of the valley operator for each eigenstate.
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, valley)

# Show bands
# save(bands, "bands.h5")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,200), colorbar=true, markercolor=:PiYG)
```