

function positionalong(lat::Lattice, v::AbstractVector; normalize=true)

    if normalize
        v = v/norm(v)
    end

    Diagonal(vec(transpose(v)*Structure.positions(lat)))
end