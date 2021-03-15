# __precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("TightBinding/types.jl") # todo : move to new module Hops
export DenseHops, SparseHops, Hops, AnyHops, AbstractHops
export hopdim, addspin, addhops!, addhops, zerokey, getzero, setzero!
export getbloch

include("TightBinding/gethops.jl") # todo : move to structure
export gethops

include("TightBinding/operators.jl")

include("TightBinding/superlattice.jl")

include("TightBinding/neighbors.jl")
export getneighborhops

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
