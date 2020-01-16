function layerprojection(lat::Lattice, d=1)
    """
    Project onto a layer
    """

    z = extrapositions(lat, "z")

    zmin = minimum(z)
    zmax = maximum(z)

    z = 2.0 .* ( (z .- zmin) ./ (zmax-zmin) .- 0.5 )

    kron(Diagonal(z), Diagonal(ones(d)))
end