module Paths

    using ...Utils

    using LinearAlgebra

    include("Paths/labeledpoints.jl")
    export LabeledPoints

    include("Paths/discretizepath.jl")

    include("Paths/path.jl")
    export DiscretePath

end