# __precompile__()

module Geometries2D

    using Base, LinearAlgebra

    using ..Structure: twist
    using ..Structure.Paths
    using ..Structure.Lattices

    # Dictionary to construct Paths in k space
    kdict_tri = LabeledPoints(
        ["γ", "κ", "μ", "κ'"],
        [[0.0;  0.0], [1/3;  2/3], [1/2;  1/2], [2/3;  1/3]],
        ["\$\\gamma\$", "\$\\kappa\$", "\$\\mu\$", "\$\\kappa'\$"],
        ["γ", "κ", "μ", "κ'", "γ", "μ"]
    )

    # Should implement this (or a similar) form in the future:
    # kdict_tri = LabeledPoints(
    #     "γ" => (label="\$\\gamma\$", coord=[0.0;  0.0]),
    #     "κ" => (label="\$\\kappa\$", coord=[1/3;  2/3]),
    #     "μ" => (label="\$\\mu\$", coord=[1/2;  1/2]),
    #     "κ'" => (label="\$\\kappa'\$", coord=[2/3;  1/3]),
    # )

    #### Define simple honeycomb systems
    A_tri = [[cos(pi/6);  -sin(pi/6)]  [cos(pi/6);  sin(pi/6)]]
    A_hex = sqrt(3) .* A_tri
    δ_hex = [1/3; 1/3]

    triangular(a::Float64=1.0) = Lattice(
        a .* A_tri,
        zeros(2,1),
        zeros(1,1);
        extradimensions=["z"],
        specialpoints=kdict_tri
    )

    # [[cos(2*pi/6);  -sin(2*pi/6)]  [cos(2*pi/6);  sin(2*pi/6)]]
    triangular_supercell(a::Float64=1.0) = Lattice(
        a .* A_tri * [[1;1] [-1;2]], #[[2;-1] [-1;2]]
        [[0.0;0.0]  [1/3;1/3]  [-1/3;2/3]],
        [zeros(1,3); 0 1 1],
        extradimensions=["z", "sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb(a::Float64=1.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  δ_hex  ],
        [ 0.0 0.0; 0 1 ];
        extradimensions=["z", "sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb_AA(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  δ_hex [0.0;0.0] δ_hex],
        [ 0 0 z z; 0 1 0 1 ];
        extradimensions=["z", "sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb_AB(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  δ_hex   -δ_hex   δ_hex],
        [ 0 0 z z; 0 1 0 1 ];
        extradimensions=["z", "sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  -δ_hex   δ_hex   -δ_hex],
        [ 0 0 z z; 0 1 0 1 ];
        extradimensions=["z", "sublattice"],
        specialpoints=kdict_tri
    )

    # honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
    #     a .* A_tri,
    #     [[0.0;0.0]   δ_hex   δ_hex  2*δ_hex],
    #     [ 0 0 z z; 0 1 0 1 ];
    #     extradimensions=["z", "sublattice"],
    #     specialpoints=kdict_tri
    # )

    function honeycomb_twisted(N::Int, a::Float64=1.0, z::Float64=3.0)
        lat = honeycomb(a)
        return twist(lat, lat, N; z=z, m=1)
    end

    function honeycomb_twisted_ABBA(N::Int, a::Float64=1.0, z::Float64=3.0)
        lat1 = honeycomb_AB(a, z)
        lat2 = honeycomb_AB(a, z)

        return twist(lat1, lat2, N; z=z, m=1)
    end

    function honeycomb_twisted_ABAB(N::Int, a::Float64=1.0, z::Float64=3.0)
        lat1 = honeycomb_AB(a, z)
        lat2 = honeycomb_BA(a, z)

        return twist(lat1, lat2, N; z=z, m=1)
    end


    function smoothdisplaceZ!(lat, δz_even=0.055, δz_odd=0.0)
        @assert latticedim(lat) == 2
        @assert abs(dot(getA(lat)[:,1], getA(lat)[:,2])/(norm(getA(lat)[:,1])*norm(getA(lat)[:,1]))) ≈ 0.5

        # Define smooth (super)lattice functions
        G1 = getB(lat)[:,1]
        G2 = getB(lat)[:,2]
        L1 = getA(lat)[:,1]
        L2 = getA(lat)[:,2]

        fodd(x::AbstractVector) =  1/(3*√3) * (sin(2*π*dot(G1,x)) + sin(2*π*dot(G2,x)) - sin(2*π*dot(G1+G2,x)))
        feven(x::AbstractVector; p0=(L1+L2)/3) = 1/3 + 2/9 * (cos(2*π*dot(G1,x-p0)) + cos(2*π*dot(G2,x-p0)) + cos(2*π*dot(G1+G2,x-p0)))

        # Modify the lattice geometry
        XY = positions(lat)
        zcoords = extrapositions(lat, "z")

        for (i_,p)=enumerate(eachcol(XY))
            zcoords[i_] += sign(zcoords[i_]) * (δz_even * feven(p) + δz_odd * fodd(p)) #-δz/2 +  # amplitude for physical problem!
        end

        setextrapositions!(lat, "z", zcoords)
    end


end
