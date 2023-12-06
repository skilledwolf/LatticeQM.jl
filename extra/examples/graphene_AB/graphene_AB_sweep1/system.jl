using LinearAlgebra

using LatticeQM

function buildsystem(;tz=0.3, V=0.3, triple=false, spinhalf=true)

    lat = Geometries.honeycomb_AB()

    if triple
        lat = Lattices.superlattice(lat, [[1,1] [2,-1]])
        Lattices.foldPC!(lat)#; shift=[1/3,-1/3,0])
    end

    H = Operators.graphene(lat; tz=tz, ℓinter=0.08, ℓintra=0.04, cellrange=3, format=:dense, mode=(spinhalf ? :spinhalf : :nospin))
    Operators.addinterlayerbias!(H, lat, V)
    
    H = DenseHops(H)

    sz = Operators.spin(lat, "sz")
    valley = Operators.valley(lat; spinhalf=spinhalf)
    # valley = Operators.valley(lat, x->x[5]; spinhalf=true)
    posZ = Operators.positionalong(lat, [0,0,1.0]; rescale=true)

    if spinhalf
        posZ = kron(posZ, Diagonal(ones(2)))
        ops = [posZ, valley, sz]
    else
        ops = [posZ, valley]
    end

    lat, H, ops
end

