__precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Base
using LinearAlgebra
using SparseArrays
using ElasticArrays, RecursiveArrayTools
using ..LatticeTools: Lattice, positions3D, get_A_3D, atom_count, has_dimension, get_positions_in

export get_operator, get_projector,  get_hops, build_H

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("../misc/paulimatrices.jl")

include("TightBinding/bloch.jl")
include("TightBinding/hamiltonian.jl")
include("TightBinding/operators.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
