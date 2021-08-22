# Getting started

## Quickstart with docker (demo)

This option requires only [docker](https://www.docker.com). All dependencies (python, julia, jupyter) will be configured automatically inside a container with the correct versions.

If you want to quickly try out the tutorial and the examples, you can start up the docker image with
```bash
$ docker-compose up --build
```

Note that the folders `extra/examples` and `extra/tutorial` are mounted as volumes. Hence, changes made to them in the container will apply also outside the container.

## Quickstart with vs-code and docker (demo and development)

Install [docker](https://www.docker.com), [vscode](https://code.visualstudio.com), as well as the plugin `Remote - container`. When you now open the git repository in vscode, it will suggest you to run the project in the container.

## Normal installation
Make sure a recent version of `julia` is installed and execute:  
```julia
# with configured SSH key:
using Pkg; Pkg.add("git@gitlab.com:skilledwolf/LatticeQM.jl.git")
# or without SSH key:
using Pkg; Pkg.add("https://gitlab.com/skilledwolf/LatticeQM.jl.git")
```

Regularly install updates with `using Pkg; Pkg.update()`.

## Example code

```@example
using Plots
using LatticeQM

# Load a lattice geometry
lat = Geometries.honeycomb()

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
plot(bands, ylabel="\$\\varepsilon/t\$", colorbar_title="valley", size=(330,200), colorbar=true, markercolor=:PiYG)
mkpath("figures") # hide
savefig("figures/graphene_bands_zeeman.svg") # hide
```

Which plots the electronic bands of graphene with a Zeeman splitting:  
![](figures/graphene_bands_zeeman.svg)

## Contributing
If you want to help develop the package, I recommend

1. Clone it onto your hard drive, i.e., navigate to `/My/Directory/` and execute `git clone URL` (replace the directory and the URL as appropriate).

2. Add the package to Julia with `using Pkg; Pkg.develop("/My/Directory/LatticeQM")`.

3. Make your modifications, test them, then commit and push.
   
4. Do not forget to regularly do `git pull` to avoid merge conflicts.

## Uninstalling
```julia
using Pkg; Pkg.remove("LatticeQM")
```