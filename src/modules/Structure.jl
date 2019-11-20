# __precompile__()

module Structure # renamed from Structure
################################################################################
################################################################################

using Base, LinearAlgebra
using RecursiveArrayTools

using HDF5

export get_kpath

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

################################################################################
################################################################################

include("Structure/type_specialpoints.jl")

include("Structure/type_lattice.jl")
include("Structure/lattice_modify.jl")
include("Structure/lattice_properties.jl")

include("Structure/type_path.jl")
include("Structure/kpath.jl")
include("Structure/grid.jl")

################################################################################
################################################################################

include("Structure/supercell.jl")
include("Structure/twist_triangular_2D.jl")

################################################################################
################################################################################
end
