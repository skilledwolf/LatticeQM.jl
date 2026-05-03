import ..Structure
import ..Structure.Lattices
import ..Structure.Lattices: Lattice

import SparseArrays: sparse
import LinearAlgebra: Diagonal, I


function s0(lat::Lattice)
    N = Lattices.countorbitals(lat)

    # 𝟙_2 ⊗ 𝟙_N (spinful identity for N orbitals)
    Diagonal(ones(2*N))
end


function Sn(lat::Lattice, n::Vector{Float64})
    N = Lattices.countorbitals(lat)
    σn = sum(n[i] .* σs[i] for i=1:3)

    # 𝟙_N ⊗ σn — spin is the fast index, matching addspin/zeeman.
    kron(sparse(1.0I, N, N), σn)
end

getsx(lat::Lattice) = Sn(lat, [1.0, 0.0, 0.0])
getsy(lat::Lattice) = Sn(lat, [0.0, 1.0, 0.0])
getsz(lat::Lattice) = Sn(lat, [0.0, 0.0, 1.0])
getsup(lat::Lattice) = 0.5 .* (s0(lat) .+ getsz(lat))
getsdown(lat::Lattice) = 0.5 .* (s0(lat) .- getsz(lat))

spin(lat::Lattice, n::Vector) = Sn(lat,n)
function spin(lat::Lattice, name::String)
    if name == "sx"
        return getsx(lat)
    elseif name == "sy"
        return getsy(lat)
    elseif name == "sz"
        return getsz(lat)
    elseif name == "sup"
        return getsup(lat)
    elseif name == "sdown"
        return getsdown(lat)
    else
        error("Name `$name` not recognized. Must be `sx`, `sy`, `sz`, `sup`, `sdown`.")
    end
end
spin(lat::Lattice, names::Vector{String}) = [spin(lat,name) for name=names]
