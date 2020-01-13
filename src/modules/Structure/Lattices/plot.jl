using Plots

@recipe function f(lat0::Lattice, colors::Union{Nothing,Vector{Float64}}=nothing; filter=(), supercell=[0:1,0:1], markercolor=:RdYlBu, clims=:auto)
    # filter could be for example filter=("layer", z->z==0.0)

    lat = deepcopy(lat0)
    repeat!(lat, supercell)

    # Plot the layers
    if filter != ()
        indices = filterindices(lat, filter...)
        atoms = positions(lat)[:,indices]
        colors = colors[indices]
    else
        atoms = positions(lat)
    end

    if colors != nothing
        # Note that colors::Vector{Float64} must be provide a color for each site
        zcolor := repeat(colors, length(supercell[1])*length(supercell[2]))
        markercolor --> markercolor
        clims --> clims
    end

    background_color_inside --> :lightgray
    seriestype := :scatter
    aspect_ratio := :equal
    markerstrokewidth := 0
    legend --> :none
    grid --> false

    xlabel --> "x"
    ylabel --> "y"

    atoms[1,:], atoms[2,:]
end