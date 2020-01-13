using Plots
using HDF5
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

function DummySave.save!(file, data::BandData) # file should be hdf5 output stream
    write(file, "bands", data.bands)
    write(file, "observable", (data.obs != nothing) ? data.obs : 0)

    DummySave.save!(file, data.path) # let's see if that works...
end

DummySave.save(data::BandData, filename::String="bands.h5") = DummySave.save_wrapper(data, filename)

@recipe function f(data::BandData, n::Integer = 1; sharpen=0.0, quantile=0.97)
    if data.obs == nothing || n == 0
        markercolor --> :black
    else
        mycolors = data.obs[:,:,n]
        max = Statistics.quantile(abs.(mycolors)[:], quantile) # maximum(abs.(data.obs[:,:,n]))

        if sharpen > 0.0
            zmin = minimum(data.obs[:,:,n])
            zmax = maximum(data.obs[:,:,n])

            mycolors = 2.0 .* ((data.obs[:,:,n].-zmin)./(zmax-zmin) .- 0.5) # scale into [-1,1]
            mycolors = tanh.(sharpen .* mycolors)
            max = 1.0
        end

        # mycolors = data.obs[:,:,n]
        zcolor      := transpose(mycolors)#transpose(mycolors)
        markercolor --> :RdYlBu #ColorGradient([:red,:limegreen,:blue])#
        clim --> (-max,max)

#         max = Statistics.quantile(abs.(mycolors)[:], 0.97)
#         clim := (-max,max)

        # max = Statistics.quantile(filter(x->x>0, data.obs[:,:,n]), 0.95)
        # min = Statistics.quantile(filter(x->x<0, data.obs[:,:,n]), 0.05)
        # max = maximum(abs.([min,max]))
        zlim --> (-max,max)

    end
    background_color_inside --> :lightgray
    ylabel --> "Energy"
    legend := :none
    seriestype  :=  :scatter
    markersize --> 1.5
    markerstrokewidth := 0
    size --> (320,300)
    xticks := (scaleticks(data.path; start=1.0, length=float(size(data.bands)[2])), data.path.ticklabels)

    transpose(data.bands)
end