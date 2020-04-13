import LatticeQM.Structure.Lattices: Lattice, allpositions, extrapositions, repeat!, sortextraposition!, filterindices, positions, countorbitals

@recipe function f(lat0::Lattice, colors::Union{Nothing,Vector{Float64},String}=nothing; filter=(), sort=false, supercell=0:0, markercolor=:RdYlBu, clims=:auto)
    # filter could be for example filter=("layer", z->z==0.0)

    if isa(colors, String) # If colors is a string it is to be understood as an extraposition
        colors = vec(extrapositions(lat0, colors))
    end

    lat = deepcopy(lat0)
    repeat!(lat, supercell)
    if colors != nothing
        colors = vcat(fill(colors, countorbitals(lat))...)
    end

    orbitalcoordinates = allpositions(lat)[1:2,:]

    perm = (:)
    if sort != false
        perm = sortextraposition!(lat, sort)
        colors = colors[perm]
        orbitalcoordinates = orbitalcoordinates[:,perm]
    end

    # Plot the layers
    if filter != ()
        indices = filterindices(lat, filter...)
        orbitalcoordinates = orbitalcoordinates[:,indices]
        colors = colors[indices]
    end

    if colors != nothing
        # Note that colors::Vector{Float64} must be provide a color for each site
        marker_z := colors
        markercolor --> markercolor
        clims --> clims
    end

    background_color_inside --> :lightgray
    seriestype := :scatter
    aspect_ratio := :equal
    markerstrokewidth := 0
    legend --> :none
    grid --> false

    xguide --> "x"
    yguide --> "y"

    orbitalcoordinates[1,:], orbitalcoordinates[2,:]
end