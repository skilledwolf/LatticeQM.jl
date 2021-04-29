import ..Structure.Lattices
import ..Structure.Lattices: Lattice

import LinearAlgebra: norm, Diagonal, I
import Statistics: mean

function positionalong(lat::Lattice, i::Integer; kwargs...)
    positionalong(lat, Lattices.basis(lat,i); kwargs...)
end


function positionalong(lat::Lattice, v::AbstractVector; rescale=false, normalize=true, center=true)

    if normalize
        v = v/norm(v)
    end

    M = Diagonal(vec(transpose(v)*Lattices.positions(lat)))

    if center
        M -= mean(diag(M))*I
    end

    if rescale
        M = 2 .*(M-minimum(M)*I)./(maximum(M)-minimum(M)) - 1.0*I
    end

    M
end