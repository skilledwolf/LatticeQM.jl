__precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Base
using LinearAlgebra
using SparseArrays
using ElasticArrays, RecursiveArrayTools

import ..Structure: Lattice, positions3D, get_A_3D, atom_count, has_dimension, get_positions_in

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export get_operator, get_projector,  get_hops, get_Hamiltonian

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("../misc/paulimatrices.jl")
include("TightBinding/types.jl")

include("TightBinding/bloch.jl")
include("TightBinding/hamiltonian.jl")
include("TightBinding/operators.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
