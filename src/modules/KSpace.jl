__precompile__()
module KSpace

using Base
# using ProgressMeter
using LinearAlgebra
using Arpack # for sparse solvers
using HDF5

# export names_to_path, names_to_coord, coord_to_points, construct_path
export PointDict, DiscretePath, eachpoint, scaled_ticks

####################################################################################

include("KSpace/struct_PathDict.jl")
include("KSpace/struct_DiscretePath.jl")

include("KSpace/samplepoints.jl")

####################################################################################
end
