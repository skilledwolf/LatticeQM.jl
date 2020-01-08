# __precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Base
using LinearAlgebra
using SparseArrays
using ElasticArrays, RecursiveArrayTools

import ..Structure: Lattice, positions, positionsND, positions3D, get_A, get_A_3D, atom_count, has_dimension, lattice_dim, assert_dimension, get_positions_in, regulargrid
import ..BlochTools: chemical_potential

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("../misc/paulimatrices.jl")
include("TightBinding/types.jl")
export DenseHops, SparseHops, Hops, AbstractHops, hopdim, extend_space, addhops!, addhops

include("TightBinding/bloch.jl")
include("TightBinding/hamiltonian.jl")
export get_hops, get_hamiltonian, get_bloch

include("TightBinding/operators.jl")
export get_operator, get_projector, @vectorwrap, add_chemicalpotential!, initial_guess, set_filling!

include("TightBinding/find_neighbors.jl")
export get_neighbors, find_common_neighbor

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
