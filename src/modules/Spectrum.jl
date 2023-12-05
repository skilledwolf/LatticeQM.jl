# __precompile__()

module Spectrum

    using LinearAlgebra
    # using KrylovKit

    include("Spectrum/eigen.jl")
    export dense

    include("Spectrum/types.jl")
    # include("Spectrum/expval.jl")

    include("Spectrum/spectrum.jl")
    # export getbands, spectrum
    include("Spectrum/chemicalpotential.jl")
    # export chemicalpotential

    include("Spectrum/bandgap.jl")
    include("Spectrum/groundstatesum.jl")

    include("Spectrum/fermisurface.jl")

    include("Spectrum/dos.jl")
    include("Spectrum/ldos.jl")
    include("Spectrum/density.jl")

    include("Spectrum/berry.jl")

end
