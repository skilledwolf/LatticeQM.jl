import PyPlot
const plt = PyPlot;

using Statistics: quantile

PiYG = plt.get_cmap("PiYG", 12)
colors = [PiYG(0.0), PiYG(0.125), (0.0,0.0,0.0,1.0), PiYG(0.875), PiYG(1.0)]
cm_PiBlYG = plt.matplotlib.colors.LinearSegmentedColormap.from_list("PiBlYG", colors, N=50)
BWR = plt.get_cmap("bwr",12)
colors = [BWR(0.0), BWR(0.125), (0.0,0.0,0.0,1.0), BWR(0.875), BWR(1.0)]
cm_BBR = plt.matplotlib.colors.LinearSegmentedColormap.from_list("BBR", colors, N=50)

function plt.plot(bands::LatticeQM.Spectrum.BandData, n::Integer=1; cmap=cm_BBR, clim=nothing, cmode=:symmetric, kwargs...)
    
    if size(bands.obs,3) == 0 || n==0
        bandcolors=:none
        vmin = nothing
        vmax = nothing
    else
        bandcolors=bands.obs[:,:,n]
        
        vmax = quantile(abs.(bandcolors)[:], 0.998)
        
        if cmode==:symmetric
            vmin = -vmax
        elseif cmode==:intensity
            vmin = 0
        elseif cmode==:minmax
            vmin = minimum(bandcolors)
            vmax = maximum(bandcolors)
        end
        
        if clim != nothing
            vmin = minimum(clim)
            vmax = maximum(clim)
        end
    end
    
    
    fig = plt.figure(;figsize=(3,2.5))
    ax = fig.add_subplot(1,1,1)

    ims = []
    for (i,b) in enumerate(eachrow(bands.bands))
        c = (bandcolors==:none) ? "k" : bandcolors[i,:]
        im = plt.scatter(bands.path.positions, b; s=1, c=c, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors="none", kwargs...)
        append!(ims,[im])
    end
    
    if bandcolors!=:none
        fig.colorbar(last(ims); ax=ax)
    end

    ticks = Structure.Paths.ticks(bands.path)
    ticklabels = Structure.Paths.tickslabels(bands.path)

    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.set_ylabel("Energy")
    
    ax
end