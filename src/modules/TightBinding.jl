# __precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Base
using LinearAlgebra
using SparseArrays
using ElasticArrays, RecursiveArrayTools

using ..Utils
import ..Algebra: σ0, σ1, σ2, σ3, σs

import ..Utils: regulargrid
import ..Structure.Lattices: Lattice, positions, allpositions, getA, countorbitals, hasdimension, latticedim, assertdimension, extrapositions

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("TightBinding/types.jl") # todo : move to new module Hops
export DenseHops, SparseHops, Hops, AnyHops, AbstractHops
export hopdim, addspin, addhops!, addhops
# export decidetype

include("TightBinding/bloch.jl") # todo : move to new module Hops
export get_bloch, getbloch

include("TightBinding/hamiltonian.jl") # todo : move to structure
export get_hops, get_hamiltonian # backwards compatibility
export gethops, gethamiltonian

include("TightBinding/operators.jl")

# Legacy definitions:
@legacymoved get_neighbors "Structure.neighbors"
@legacymoved find_common_neighbor "Structure.commonneighbor"
@legacymoved initial_guess "Meanfield.initialguess"
@legacymoved set_filling! "Operatorssetfilling!"
@legacymoved add_chemicalpotential! "Operatorsaddchemicalpotential!"

@legacymoved get_operator "Operators.getoperator"
@legacymoved get_projector "Operators.get_projector"
export get_operator, get_projector

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
