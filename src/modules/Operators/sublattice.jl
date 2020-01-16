function sublatticeprojection(lat::Lattice, n::Number, d=1)
    """
    typically n=0 for sublattice A and n=1 for sublattice B and d=2 for spin-1/2
    """

    sublattice = extrapositions(lat, "sublattice")

    filtered_sublattice = Array{Float64}(sublattice .== n)


    kron(Diagonal(filtered_sublattice), Diagonal(ones(d)))
end

function sublatticeprojection(lat::Lattice, n::String, args...)
    if n=="A"
        return sublattice(lat,0, args...)
    elseif n=="B"
        return sublattice(lat,1,args...)
    else
        error("Did not recognize sublattice `$n`. Must be `A` or `B` or integer.")
    end
end

