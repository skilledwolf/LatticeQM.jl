
"""
    LatticeQM 

Library for tight-binding models defined on (periodic) lattices, providing 
convenient functions to build the operators and to obtain bands, expectation values,
topological indices, linear response coefficients and mean-field solutions.

## Submodules
- Structure
- TightBinding 
- Spectrum
- Operators
- LinearResponse
- Meanfield

Any of these modules can be further explored, e.g., with `?TightBinding`

## Usage examples
See folder `examples` of the package.

"""
module LatticeQM

# using ElasticArrays
using SparseArrays, LinearAlgebra

# using KPM # this is one of my packages, removed dependency

include("modules/Utils/Utils.jl")
import .Utils
import .Utils: dense

include("modules/Eigen/Eigen.jl")
import .Eigen

### Base modules
include("modules/Structure/Structure.jl")
import .Structure
import .Structure: kpath, Geometries, Lattices
export Structure, Geometries, Lattices, kpath

include("modules/Spectrum/Spectrum.jl")
import .Spectrum
import .Spectrum: getbands
export Spectrum, getbands

include("modules/TightBinding/TightBinding.jl")
import .TightBinding
import .TightBinding: DenseHops, SparseHops, Hops, hopdim, addhops!, addhops
export TightBinding
export DenseHops, SparseHops, Hops, hopdim, addhops!, addhops
export dense, sparse

include("modules/Floquet/Floquet.jl")
using .Floquet
export Floquet

include("modules/Operators/Operators.jl")
import .Operators
import .Operators: gethops, getoperator, setfilling! #, getprojector
export Operators
export gethops, getoperator, setfilling! #, getprojector

include("modules/LinearResponse/LinearResponse.jl")
import .LinearResponse
export LinearResponse

include("modules/Meanfield/Meanfield.jl")
import .Meanfield
import .Meanfield: solveselfconsistent, solvehartreefock, HartreeFock, initialguess
export Meanfield
export solveselfconsistent, solvehartreefock, HartreeFock, initialguess

include("modules/Superconductivity/Superconductivity.jl")
using .Superconductivity
export Superconductivity
export BdGOperator

include("modules/Plotting/Plotting.jl")
using .Plotting

end # module LatticeQM
