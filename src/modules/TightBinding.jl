# __precompile__()

module TightBinding
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include("TightBinding/types.jl") # todo : move to new module Hops
export DenseHops, SparseHops, Hops, AbstractHops
export hopdim, addspin, addhops!, addhops, zerokey, getzero, setzero!
export dense, sparse


include("TightBinding/operators.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
end
