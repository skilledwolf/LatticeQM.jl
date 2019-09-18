
using ..Structure: DiscretePath, scaled_ticks

export bandplot, bandplot!

function bandplot(bands::Matrix{Float64}, ks::DiscretePath; μ=0.0, zcolor=nothing, markersize=3, cmap=:RdBu, ylims=nothing, plotsize=(300,480), kwargs...)

    p=plot()
    bandplot!(p, bands, ks; μ=μ, zcolor=zcolor, markersize=markersize, cmap=cmap, ylims=ylims, plotsize=plotsize, kwargs...)

    p
end

function bandplot!(p, bands::Matrix{Float64}, ks::DiscretePath; μ=0.0, zcolor=nothing, markersize=3, cmap=:RdBu, ylims=nothing, plotsize=(300,480), kwargs...)

    if zcolor == nothing
        scatter!(p,transpose(bands),
            legend=:none,
            m=(:black, 0.8, Plots.stroke(0, :none)),
            markersize=markersize,
            size=plotsize,
            ylims=ylims,
            # marker=".",
            ylabel="\$\\varepsilon/t\$",
            kwargs...
        )
    else
        scatter!(p,transpose(bands),
            zcolor=zcolor,
            legend=:none,
            m=(:RdBu, 0.8, Plots.stroke(0, :none)),
            markersize=markersize,
            size=plotsize,
            ylims=ylims,
            # marker=".",
            ylabel="\$\\varepsilon/t\$",
            kwargs...
        )
    end

    xticks!(p, scaled_ticks(ks; start=1.0, length=float(size(bands)[2])), ks.ticklabels)

    nothing
end
