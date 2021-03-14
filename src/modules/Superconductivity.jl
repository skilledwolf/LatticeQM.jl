
module Superconductivity

    using LinearAlgebra

    include("Superconductivity/types.jl")
    export BdGOperator

    include("Superconductivity/operators.jl")

    include("Superconductivity/green.jl")

    include("Superconductivity/hartreefock.jl")

end