# __precompile__()
module Materials # proposal: rename this module to LatticeOperators

using LinearAlgebra
using SparseArrays

using ..Utils
using ..Utils: @scalar2vector
using ..Algebra: σ0, σ1, σ2, σ3, σs

using ..Structure: regulargrid
using ..Structure.Lattices: Lattice, assertdimension, hasdimension, extrapositions, countatoms, positions, getA, allpositions, latticedim
using ..TightBinding
using ..Spectrum: chemicalpotential

include("Materials/nearestneighbor.jl")

include("Materials/chemicalpotential.jl")
export add_chemicalpotential!, initial_guess, set_filling!

include("Materials/zeeman.jl")
export getzeeman, addzeeman!

include("Materials/graphene.jl")

####################################################################################
end
