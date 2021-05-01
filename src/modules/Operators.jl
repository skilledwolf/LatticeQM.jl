# __precompile__()
module Operators # proposal: rename this module to LatticeOperators

using LinearAlgebra
using SparseArrays

using ..Utils
using ..Utils: @scalar2vector

using ..Structure
using ..Structure.Lattices
using ..Structure: regulargrid
using ..TightBinding
using ..Spectrum: chemicalpotential

include("Operators/gethops.jl")
export gethops

include("Operators/superlattice.jl")

include("Operators/getoperator.jl")
export getoperator, getprojector

include("Operators/trace.jl")
export trace, expval, magnetization

include("Operators/spin.jl")

include("Operators/neighbors.jl")

include("Operators/nearestneighbor.jl")

include("Operators/chemicalpotential.jl")
export addchemicalpotential!, setfilling!

include("Operators/densitymatrix.jl")

include("Operators/zeeman.jl")

include("Operators/current.jl")

include("Operators/peierls.jl")

include("Operators/interactions.jl")
# export gethubbard, getcappedyukawa

include("Operators/graphene.jl") # adds 1.5 seconds to the load time due to PyCall

include("Operators/haldanelike.jl")

include("Operators/sublattice.jl")
include("Operators/layer.jl")

include("Operators/position.jl")

####################################################################################
end
