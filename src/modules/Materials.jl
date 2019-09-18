# __precompile__()
module Materials

using LinearAlgebra

using ..Structure: Lattice, assert_dimension, get_positions_in, atom_count, positions, lattice_dim
using ..TightBinding: get_hamiltonian, get_bloch, get_hops, extend_space

include("../misc/paulimatrices.jl")

include("Materials/graphene.jl")



####################################################################################
end
