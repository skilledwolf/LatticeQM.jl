import Statistics: quantile
import LatticeQM.Spectrum: BandData
import LatticeQM.Structure.Paths: scaleticks

@recipe function f(data::BandData, n::Integer = 1; sharpen=0.0, quant=0.97)
    if data.obs == nothing || n == 0
        markercolor --> :black
    else
        mycolors = data.obs[:,:,n]
        max = quantile(abs.(mycolors)[:], quant) # maximum(abs.(data.obs[:,:,n]))

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
    xticks --> (scaleticks(data.path; start=1.0, length=float(size(data.bands)[2])), data.path.ticklabels)

    transpose(data.bands)
end