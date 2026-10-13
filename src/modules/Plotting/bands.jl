import Statistics: quantile
import LatticeQM.Spectrum
import LatticeQM.Structure.Paths

@recipe function f(data::Spectrum.BandData, n::Integer=1; sharpen=0.0, climits=nothing, csymmetric=true, cquantile=0.97)
    bands = data.bands
    if data.obs == nothing || n == 0 || size(data.obs,3) == 0
        markercolor --> :black
    else
        mycolors = data.obs[:,:,n]
        # Draw colors in order without changing the stored band/observable pairs.
        perm = sortperm(mycolors; dims=1)
        bands = bands[perm]
        mycolors = mycolors[perm]

        max = quantile(abs.(mycolors)[:], cquantile) # maximum(abs.(data.obs[:,:,n]))
        # max = maximum(abs.(data.obs[:,:,n]))
        min = (csymmetric ? -max : quantile(abs.(mycolors)[:], 1-cquantile))


        if sharpen > 0.0
            zmin, zmax = extrema(mycolors)

            mycolors = zmax == zmin ? zero(mycolors) :
                tanh.(sharpen .* (2 .* (mycolors .- zmin) ./ (zmax-zmin) .- 1))
            max = 1.0
            min = -1.0
        end
        # mycolors = data.obs[:,:,n]
        marker_z      := transpose(mycolors)#transpose(mycolors)
        markercolor --> :diverging_bkr_55_10_c35_n256 #:RdYlBu #ColorGradient([:red,:limegreen,:blue])#

        clims --> ((climits==nothing) ? (min-sqrt(eps()), max+sqrt(eps())) : climits)

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
    xticks --> (Paths.scaleticks(data.path; start=1.0, length=float(size(data.bands)[2])), Paths.ticklabels(data.path))

    transpose(bands)
end
