# __precompile__()
module Operators # proposal: rename this module to LatticeOperators

using LinearAlgebra
using SparseArrays

using ..Utils
using ..Utils: @scalar2vector
using ..Algebra

using ..Structure
using ..Utils: regulargrid
using ..Structure.Lattices
using ..TightBinding
using ..Spectrum: chemicalpotential

include("Operators/getoperator.jl")
export getoperator, getprojector

include("Operators/trace.jl")
export trace, expval, magnetization

include("Operators/spin.jl")

include("Operators/nearestneighbor.jl")

include("Operators/chemicalpotential.jl")
export addchemicalpotential!, setfilling!

include("Operators/zeeman.jl")
export getzeeman, addzeeman!

include("Operators/current.jl")

include("Operators/graphene.jl")

include("Operators/sublattice.jl")
include("Operators/layer.jl")

####################################################################################
end
