import ..Structure
import ..Structure.Lattices: Lattice

import SparseArrays: spzeros
import LinearAlgebra: Diagonal


function s0(lat::Lattice)
    N = Structure.countorbitals(lat)

    # d.σ ⊗ 𝟙_N
    Diagonal(ones(2*N))
end

@legacyalias Sn S_n
function Sn(lat::Lattice, n::Vector{Float64})
    N = Structure.countorbitals(lat)

    # d.σ ⊗ 𝟙_N
    mat = spzeros(ComplexF64, 2*N, 2*N)
    σn = sum(n[i] .* σs[i] for i=1:3)

    @simd for i = 1:2:2*N
        mat[i:i+1, i:i+1] .= σn
    end

    mat
end

getsx(lat::Lattice) = Sn(lat, [1.0, 0.0, 0.0])
getsy(lat::Lattice) = Sn(lat, [0.0, 1.0, 0.0])
getsz(lat::Lattice) = Sn(lat, [0.0, 0.0, 1.0])
getsup(lat::Lattice) = 0.5 .* (S0(lat) .+ getsz(lat))
getsdown(lat::Lattice) = 0.5 .* (S0(lat) .- getsz(lat))

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
