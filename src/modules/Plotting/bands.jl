
using ..KSpace: DiscretePath, scaled_ticks

export bandplot, bandplot!

function bandplot(bands::Matrix{Float64}, ks::DiscretePath; ylims=nothing, plotsize=(300,480))

    p=plot()
    bandplot!(p, bands, ks; ylims=ylims, plotsize=plotsize)

    p
end

function bandplot!(p, bands::Matrix{Float64}, ks::DiscretePath; ylims=nothing, plotsize=(300,480))

    plot!(p,transpose(bands),
        legend=:none,
        c=:cornflowerblue,
        size=plotsize,
        ylims=ylims,
        ylabel="\$\\varepsilon/t\$"
    )
    xticks!(p, scaled_ticks(ks; start=1.0, length=float(size(bands)[2])), ks.ticklabels)

    nothing
end
