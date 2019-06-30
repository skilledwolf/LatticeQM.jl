
using ..KSpace: DiscretePath, scaled_ticks

export bandplot, bandplot!

function bandplot(bands::Matrix{Float64}, ks::DiscretePath)

    p=plot()
    bandplot!(p, bands, ks)

    p
end

function bandplot!(p, bands::Matrix{Float64}, ks::DiscretePath)

    plot!(p,transpose(bands),
        legend=:none,
        c=:cornflowerblue,
        size=(300,480),
        ylims=(-0.3,0.3),
        ylabel="\$\\varepsilon/t\$"
    )
    xticks!(p, scaled_ticks(ks; start=1.0, length=float(size(bands)[2])), ks.ticklabels)

    nothing
end
