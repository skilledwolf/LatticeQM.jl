# __precompile__()
module Materials

using LinearAlgebra

using ..Structure: Lattice, assert_dimension, has_dimension, get_positions_in, atom_count, positions, get_A, positions3D, get_A_3D, positionsND, lattice_dim
using ..TightBinding: get_hamiltonian, get_bloch, get_hops, extend_space, get_neighbors, find_neighbors, find_common_neighbor, add_hoppings!

include("../misc/paulimatrices.jl")

include("Materials/generic.jl")
include("Materials/graphene.jl")

####################################################################################
end
