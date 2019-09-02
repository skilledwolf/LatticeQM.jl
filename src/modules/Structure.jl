__precompile__()

module Structure # renamed from Structure
################################################################################
################################################################################

using Base, LinearAlgebra
using RecursiveArrayTools

using HDF5

################################################################################
################################################################################

include("Structure/specialpoints.jl")
include("Structure/kpath.jl")
include("Structure/grid.jl")

include("Structure/lattice.jl")
include("Structure/lattice_modify.jl")
include("Structure/lattice_properties.jl")

include("Structure/export.jl")

################################################################################
################################################################################

include("Structure/supercell.jl")
include("Structure/trianguar_twist_2D.jl")

################################################################################
################################################################################
end
