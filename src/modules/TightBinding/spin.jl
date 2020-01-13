function S0(lat::Lattice)
    N = countatoms(lat)

    # d.σ ⊗ 𝟙_N
    Diagonal(ones(2*N))
end

@legacyalias Sn S_n
function Sn(lat::Lattice, n::Vector{Float64})
    N = countatoms(lat)

    # d.σ ⊗ 𝟙_N
    mat = spzeros(Complex, 2*N, 2*N)
    σn = sum(n[i] .* σs[i] for i=1:3)

    @simd for i = 1:2:2*N
        mat[i:i+1, i:i+1] .= σn
    end

    mat
end

SX(lat::Lattice) = S_n(lat, [1.0, 0.0, 0.0])
SY(lat::Lattice) = S_n(lat, [0.0, 1.0, 0.0])
SZ(lat::Lattice) = S_n(lat, [0.0, 0.0, 1.0])

Sup(lat::Lattice) = 0.5 .* (S0(lat) .+ SZ(lat))
Sdown(lat::Lattice) = 0.5 .* (S0(lat) .- SZ(lat))

MX = SX
MY = SY
MZ = SZ