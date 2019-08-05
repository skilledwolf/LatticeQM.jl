__precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Base
using LinearAlgebra
using SparseArrays
using ElasticArrays, RecursiveArrayTools
using ..LatticeTools: Lattice, positions3D, get_A_3D, atom_count, has_dimension

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("../misc/paulimatrices.jl")

include("TightBinding/build_BlochH.jl")
include("TightBinding/build_H.jl")
include("TightBinding/build_operators.jl")
include("TightBinding/build_meanfield_op.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
