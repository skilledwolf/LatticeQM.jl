using Plots
using HDF5
using ..Structure: scaled_ticks

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

function DummySave.save!(file, data::BandData) # file should be hdf5 output stream
    write(file, "bands", data.bands)
    write(file, "observable", (data.obs != nothing) ? data.obs : 0)

    DummySave.save!(file, data.path) # let's see if that works...
end

DummySave.save(data::BandData, filename::String="bands.h5") = DummySave.save_wrapper(data, filename)

@recipe function f(data::BandData, n::Integer = 1)
    if data.obs == nothing || n == 0
        markercolor --> :blue
    else
        zcolor      := transpose(data.obs[:,:,n])
        markercolor := :RdYlBu
    end
    background_color_inside --> :lightgray
    ylabel --> "Energy"
    legend := :none
    seriestype  :=  :scatter
    markersize --> 1.5
    markerstrokewidth := 0
    size --> (320,300)
    xticks := (scaled_ticks(data.path; start=1.0, length=float(size(data.bands)[2])), data.path.ticklabels)

    transpose(data.bands)
end

################################################################################
################################################################################
