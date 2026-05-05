# __precompile__()
module Operators # proposal: rename this module to LatticeOperators

using LinearAlgebra
using SparseArrays

using ..Utils
using ..Utils: @scalar2vector

import ..Structure
import ..Structure.Lattices
import ..Structure: regulargrid
import ..TightBinding
import ..TightBinding: Hops
import ..Spectrum: chemicalpotential

include("gethops.jl")
export gethops

include("superlattice.jl")

include("getoperator.jl")
export getoperator, getprojector

include("trace.jl")
export trace, expval, magnetization

include("spin.jl")

include("neighbors.jl")

include("nearestneighbor.jl")

include("chemicalpotential.jl")
export addchemicalpotential!, setfilling!

include("densitymatrix.jl")

include("zeeman.jl")

include("current.jl")

include("peierls.jl")

include("interactions.jl")
# export gethubbard, getcappedyukawa

include("graphene.jl")

include("haldanelike.jl")

include("sublattice.jl")
include("layer.jl")

include("position.jl")

####################################################################################
end
