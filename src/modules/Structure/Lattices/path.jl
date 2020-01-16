
function Paths.DiscretePath(lat::Lattice; kwargs...)
    named_path = lat.specialpoints.defaultpath

    Paths.DiscretePath(lat, named_path; kwargs...)
end

function Paths.DiscretePath(lat::Lattice, named_path::Vector{String}; kwargs...)
    Paths.DiscretePath(lat.specialpoints, named_path; B=getB(lat), kwargs...)
end

function Paths.DiscretePath(lat::Lattice, kdict::LabeledPoints, named_path::Vector{String}; kwargs...)
    Paths.DiscretePath(kdict, named_path; B=getB(lat), kwargs...)
end

@legacyalias kpath get_kpath
function kpath(lat::Lattice, args...; kwargs...)
    """
    alias of DiscretePath(lat).
    """

    Paths.DiscretePath(lat, args...; kwargs...)
end