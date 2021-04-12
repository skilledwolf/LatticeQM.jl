function reduceto0D(lat::Lattice, N::Int)
    d=latticedim(lat)
    v = ones(Int,d)
    v[d]=N
    reduceto0D(lat, v)
end
reduceto0D(lat::Lattice, v::AbstractVector{Int}) = reduceto0D(lat, Diagonal(v))
function reduceto0D(lat::Lattice, M::AbstractMatrix{Int})
    lat1 = superlattice(lat, M)
    lat1.latticedim = 0
    kdict = Paths.LabeledPoints(
        ["Γ"],
        [[0.0]],
        ["\$\\Gamma\$"],
        ["Γ"]
    )
    lat1.specialpoints=kdict

    lat1
end


function reduceto1D(lat::Lattice, N::Int)
    d = latticedim(lat)
    v = ones(Int,d)
    v[d] = N
    reduceto1D(lat, v)
end
reduceto1D(lat::Lattice, v::AbstractVector{Int}) = reduceto1D(lat, Matrix(Diagonal(v)))
function reduceto1D(lat::Lattice, M::AbstractMatrix{Int})
    lat1 = superlattice(lat, M)
    lat1.latticedim = 1
    kdict = Paths.LabeledPoints(
        ["X1", "Γ", "X2", "Γ2"],
        [[-0.5], [0.0], [0.5], [1.0]],
        ["\$-G/2\$", "\$0\$", "\$G/2\$", "\$G\$"],
        ["X1", "Γ", "X2"]
    )
    lat1.specialpoints=kdict

    lat1
end