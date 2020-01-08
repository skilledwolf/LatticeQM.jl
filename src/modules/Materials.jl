# __precompile__()
module Materials

using LinearAlgebra

using ..Structure: Lattice, assert_dimension, has_dimension, get_positions_in, atom_count, positions, get_A, positions3D, get_A_3D, positionsND, lattice_dim, regulargrid
using ..TightBinding
using ..BlochTools: chemical_potential

include("../misc/paulimatrices.jl")

include("Materials/generic.jl")
include("Materials/graphene.jl")

####################################################################################
end
