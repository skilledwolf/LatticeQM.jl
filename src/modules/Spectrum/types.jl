# using Plots
# using HDF5
using DelimitedFiles
import ..DummySave

import ..Structure.Paths: DiscretePath

# This file is meant to contain "result types", i.e.,
# containers for data, methods to store them and to visualize them.

################################################################################
################################################################################

"""
    BandData

Struct to store BandData obtained from method `getbands(...)`.
The fields `BandData.bands` is a matrix of eigenvalues at each point along a path,
the field `BandData.obs` is a matrix of expectation values and `BandData.path` is
a DiscretePath object (contains discrete points and point labels).

Can be saved conveniently with `savedlm(bands; path="data")` and plotted with `plot(bands)`.``

"""
mutable struct BandData
    bands::Matrix{Float64}
    obs::Union{Nothing, Array{Float64,3}}
    path::DiscretePath
end

function Base.show(io::IO, bands::BandData)
    println(io, "Number of bands:      ", size(bands.bands,1))
    println(io, "Number of k-points:   ", size(bands.bands,2))
    
    if isa(bands.obs,Nothing)
        println(io, "No observables.")
    else
        println(io, "Number of observables: ", size(bands.obs,3))
    end
    show(io, bands.path)
end

using ..Structure.Paths: ticks, ticklabels, scaleticks

function DummySave.savedlm(bands::BandData; path="data", suffix="")
    mkpath(path)
    writedlm(joinpath(path, "kticks$suffix.out"), scaleticks(bands.path; start=1.0, length=float(size(bands.bands)[2])))
    writedlm(joinpath(path, "kticklabels$suffix.out"), ticklabels(bands.path))
    writedlm(joinpath(path, "bands$suffix.out"), bands.bands)
    if bands.obs != nothing
        for k=1:size(bands.obs, 3)
            writedlm(joinpath(path, "bandcolors$suffix$k.out"), bands.obs[:,:,k])
        end
    end
end

function DummySave.save!(file, data::BandData) # file should be hdf5 output stream
    write(file, "bands", data.bands)
    write(file, "observable", (data.obs != nothing) ? data.obs : 0)

    DummySave.save!(file, data.path) # let's see if that works...
end

DummySave.save(data::BandData, filename::String="bands.h5") = DummySave.save_wrapper(data, filename)
