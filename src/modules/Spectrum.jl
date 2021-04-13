# __precompile__()

module Spectrum

    export getbands, spectrum, chemicalpotential

    include("Spectrum/types.jl")

    include("Spectrum/expval.jl")
    include("Spectrum/spectrum.jl")
    include("Spectrum/chemicalpotential.jl")
    include("Spectrum/bandgap.jl")
    include("Spectrum/groundstatesum.jl")

    include("Spectrum/dos.jl")
    include("Spectrum/ldos.jl")

    include("Spectrum/berry.jl")

    #include("BlochTools/meanfield.jl")
end
