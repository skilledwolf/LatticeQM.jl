import ..Structure
import ..Structure.Lattices: Lattice

import LinearAlgebra: norm, Diagonal


function positionalong(lat::Lattice, v::AbstractVector; normalize=true)

    if normalize
        v = v/norm(v)
    end

    Diagonal(vec(transpose(v)*Structure.positions(lat)))
end