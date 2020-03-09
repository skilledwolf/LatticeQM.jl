using Plots

@recipe function f(lat0::Lattice, colors::Union{Nothing,Vector{Float64},String}=nothing; filter=(), sort=false, supercell=[0:1,0:1], markercolor=:RdYlBu, clims=:auto)
    # filter could be for example filter=("layer", z->z==0.0)

    if isa(colors, String) # If colors is a string it is to be understood as an extraposition
        colors = extrapositions(lat0, colors)
    end

    lat = deepcopy(lat0)
    repeat!(lat, supercell)

    perm = (:)
    if sort != false
        perm = sortextraposition!(lat, sort)
    end

    # Plot the layers
    if filter != ()
        indices = filterindices(lat, filter...)
        orbitalcoordinates = positions(lat)[:,indices]
        colors = colors[indices]
    else
        orbitalcoordinates = positions(lat)
    end

    if colors != nothing
        # Note that colors::Vector{Float64} must be provide a color for each site
        zcolor := vcat(fill(colors, countorbitals(lat0))...)[perm]
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

    orbitalcoordinates[1,:], orbitalcoordinates[2,:]
end