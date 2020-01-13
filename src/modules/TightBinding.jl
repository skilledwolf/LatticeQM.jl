# __precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Base
using LinearAlgebra
using SparseArrays
using ElasticArrays, RecursiveArrayTools

using ..Utils
import ..Algebra: σ0, σ1, σ2, σ3, σs

import ..Structure: regulargrid
import ..Structure.Lattices: Lattice, positions, allpositions, getA, countatoms, hasdimension, latticedim, assertdimension, extrapositions
import ..Spectrum: chemicalpotential

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("TightBinding/types.jl") # todo : move to new module Hops
export DenseHops, SparseHops, Hops, AnyHops, AbstractHops
export hopdim, addspin, addhops!, addhops

include("TightBinding/bloch.jl") # todo : move to new module Hops
export get_bloch, getbloch

include("TightBinding/hamiltonian.jl") # todo : move to structure
export get_hops, get_hamiltonian # backwards compatibility
export gethops, gethamiltonian

include("TightBinding/operators.jl")
export getoperator, getprojector

include("TightBinding/spin.jl")

@legacymoved get_neighbors "Structure.neighbors"
@legacymoved find_common_neighbor "Structure.commonneighbor"
export get_neighbors, find_common_neighbor

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
