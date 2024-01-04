
module Superconductivity

    using LinearAlgebra

    include("types.jl")
    export BdGOperator

    include("operators.jl")

    include("green.jl")

    include("hartreefock.jl")

end