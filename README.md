# LatticeQM.jl

This is a development repository for julia package that deals with lattice structures and tight-binding models.

Not intended for public use yet. Much of the documentation is still missing.  
To get started, check out the notebooks in the folder `tutorial`.  
If you want to know more, get in touch with me.


## Installation
```julia
using Pkg; Pkg.add("https://github.com/skilledwolf/KPM.jl.git")
using Pkg; Pkg.add("https://github.com/skilledwolf/LatticeQM.jl.git")
```
or for the SSH version (needs SSH key to be set up)
```julia
using Pkg; Pkg.add("git@github.com:skilledwolf/KPM.jl.git")
using Pkg; Pkg.add("git@github.com:skilledwolf/LatticeQM.jl.git")
```

Regurlarly install updates with `using Pkg; Pkg.update()`.

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
valley = Operators.valley(lat; spinhalf=true)

# Get bands along (default) high-symmetry path
# and compute the expectation value of the valley operator for each eigenstate.
ks = kpath(lat; num_points=200)
bands = getbands(hops, ks, valley)

# Show bands
# save(bands, "bands.h5")
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,200), colorbar=true, markercolor=:PiYG)
```

## Contributing
If you want to help develop the package, I recommend

1. Clone it onto your hard drive, i.e., navigate to `/My/Directory/` and execute `git clone URL` (replace the directory and the URL as appropriate).

2. Add the package to Julia with `using Pkg; Pkg.dev("/My/Directory/LatticeQM")`.

3. Make your modifications, test them, then commit and push.
   
4. Do not forget to regularly do `git pull` to avoid merge conflicts.

## Uninstalling
```julia
using Pkg; Pkg.remove("LatticeQM")
```
