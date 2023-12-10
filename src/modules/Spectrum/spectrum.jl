# __precompile__()

module Spectrum

    using LinearAlgebra
    # using KrylovKit

    import LatticeQM.Eigen

    include("types.jl")
    # include("expval.jl")

    include("bands.jl")
    # export getbands, spectrum
    include("chemicalpotential.jl")
    # export chemicalpotential

    include("bandgap.jl")
    include("groundstatesum.jl")

    include("fermisurface.jl")

    include("dos.jl")
    include("ldos.jl")
    include("density.jl")

    include("berry.jl")

end
