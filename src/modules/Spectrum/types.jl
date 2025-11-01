
import ..Structure.Paths: DiscretePath

# This file is meant to contain "result types", i.e.,
# containers for data, methods to store them and to visualize them.

################################################################################
################################################################################

import SharedArrays

"""
    BandData

Struct to store BandData obtained from method `getbands(...)`.
The fields `BandData.bands` is a matrix of eigenvalues at each point along a path,
the field `BandData.obs` is a matrix of expectation values and `BandData.path` is
a DiscretePath object (contains discrete points and point labels).

Can be saved conveniently with `savedlm(bands; path="data")` and plotted with `plot(bands)`.``

"""
mutable struct BandData{T} # T = DiscretePath in practice
    bands::Matrix{Float64}
    obs::Array{Float64,3}
    path::T

    BandData(bands::AbstractMatrix, obs::AbstractArray, path::T) where {T} = new{T}(SharedArrays.sdata(bands), SharedArrays.sdata(obs), path)
end

function Base.show(io::IO, bands::BandData)
    println(io, "Number of bands:      ", size(bands.bands,1))
    println(io, "Number of k-points:   ", size(bands.bands,2))
    
    if size(bands.obs,3)==0
        println(io, "No observables.")
    else
        println(io, "Number of observables: ", size(bands.obs,3))
    end
    show(io, bands.path)
end

using ..Structure.Paths: ticks, ticklabels, scaleticks


using DelimitedFiles

DelimitedFiles.writedlm(bands::BandData) = DelimitedFiles.writedlm(".", bands)

function DelimitedFiles.writedlm(path::String, bands::BandData; suffix="")
    @assert ispath(path) "Error: '$path' is not a valid path."
    writedlm(joinpath(path, "kticks$suffix.out"), scaleticks(bands.path; start=1.0, length=float(size(bands.bands)[2])))
    writedlm(joinpath(path, "kticklabels$suffix.out"), ticklabels(bands.path))
    writedlm(joinpath(path, "bands$suffix.out"), bands.bands)
    if size(bands.obs, 3)>0
        for k=1:size(bands.obs, 3)
            writedlm(joinpath(path, "bandcolors$suffix$k.out"), bands.obs[:,:,k])
        end
    end
end
