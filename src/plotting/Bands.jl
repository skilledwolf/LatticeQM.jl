import Statistics: quantile
import LatticeQM.Spectrum: BandData
import LatticeQM.Structure.Paths: scaleticks

@recipe function f(data::BandData, n::Integer = 1; sharpen=0.0, clims=nothing, csymmetric=true, cquantile=0.97)
    if data.obs == nothing || n == 0
        markercolor --> :black
    else
        mycolors = data.obs[:,:,n]
        max = quantile(abs.(mycolors)[:], cquantile) # maximum(abs.(data.obs[:,:,n]))
        # max = maximum(abs.(data.obs[:,:,n]))

        if csymmetric
            min = -max
        else
            min = quantile(abs.(mycolors)[:], 1-cquantile)
        end

        if sharpen > 0.0
            zmin = minimum(data.obs[:,:,n])
            zmax = maximum(data.obs[:,:,n])

            mycolors = 2.0 .* ((data.obs[:,:,n].-zmin)./(zmax-zmin) .- 0.5) # scale into [-1,1]
            mycolors = tanh.(sharpen .* mycolors)
            max = 1.0
        end
        # mycolors = data.obs[:,:,n]
        marker_z      := transpose(mycolors)#transpose(mycolors)
        markercolor --> :diverging_bkr_55_10_c35_n256 #:RdYlBu #ColorGradient([:red,:limegreen,:blue])#

        if clims==nothing
            clims --> (x->(min-sqrt(eps()),max+sqrt(eps())))
        else
            clims --> (x->clims)
        end

#         max = Statistics.quantile(abs.(mycolors)[:], 0.97)
#         clims := (-max,max)

        # max = Statistics.quantile(filter(x->x>0, data.obs[:,:,n]), 0.95)
        # min = Statistics.quantile(filter(x->x<0, data.obs[:,:,n]), 0.05)
        # max = maximum(abs.([min,max]))
        zlims --> (min,max)

    end
    background_color_inside --> :lightgray
    yguide --> "Energy"
    legend := :none
    seriestype  :=  :scatter
    markersize --> 1.5
    markerstrokewidth := 0
    size --> (320,300)
    xticks --> (scaleticks(data.path; start=1.0, length=float(size(data.bands)[2])), data.path.ticklabels)

    transpose(data.bands)
end