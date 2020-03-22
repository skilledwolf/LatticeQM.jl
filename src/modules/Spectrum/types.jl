# using Plots
# using HDF5
using DelimitedFiles
using Statistics
using ..Structure.Paths: scaleticks

using ..DummySave

# This file is meant to contain "result types", i.e.,
# containers for data, methods to store them and to visualize them.

################################################################################
################################################################################

mutable struct BandData
    bands::Matrix{Float64}
    obs::Union{Nothing, Array{Float64,3}}
    path::DiscretePath
end

function DummySave.savedlm(bands::BandData; path="data", suffix="")
    mkpath(path)
    writedlm(joinpath(path, "kticks$suffix.out"), scaleticks(bands.path; start=1.0, length=float(size(bands.bands)[2])))
    writedlm(joinpath(path, "kticklabels$suffix.out"), bands.path.ticklabels)
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
