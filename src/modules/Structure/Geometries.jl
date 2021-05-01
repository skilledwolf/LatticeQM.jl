# __precompile__()


"""
    Geometries

Provides predefined lattice objects (such as two-dimensional honeycomb lattice).

### Example
```julia
import Structure.Geometries

lat = Geometries.honeycomb_twisted(11)
plot(lat, 3; supercell=[0:1,0:1])
```
"""
module Geometries

include("Geometries/2d.jl")

end