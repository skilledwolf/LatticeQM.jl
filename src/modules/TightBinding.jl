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

export get_operator, get_projector, get_hops, get_hamiltonian, get_bloch, get_dense, hopdim, initial_guess, set_filling!

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("../misc/paulimatrices.jl")
include("TightBinding/types.jl")

include("TightBinding/bloch.jl")
include("TightBinding/hamiltonian.jl")
include("TightBinding/operators.jl")

include("TightBinding/find_neighbors.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
