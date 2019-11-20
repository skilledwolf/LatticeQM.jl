__precompile__()

module Geometries2D

    using Base, LinearAlgebra

    import ..Structure: Lattice, PointDict, has_dimension, get_positions_in, atom_count, positions, twist_triangular_2D

    # Dictionary to construct Paths in k space
    kdict_tri = PointDict(
        ["γ", "κ", "μ", "κ'"],
        [[0.0;  0.0], [1/3;  2/3], [1/2;  1/2], [2/3;  1/3]],
        ["\$\\gamma\$", "\$\\kappa\$", "\$\\mu\$", "\$\\kappa'\$"],
        ["γ", "κ", "μ", "κ'", "γ"]
    )

    # Should implement this (or a similar) form in the future:
    # kdict_tri = PointDict(
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
        highsymmetrypoints=kdict_tri
    )

    # [[cos(2*pi/6);  -sin(2*pi/6)]  [cos(2*pi/6);  sin(2*pi/6)]]
    triangular_supercell(a::Float64=1.0) = Lattice(
        a .* A_tri * [[1;1] [-1;2]], #[[2;-1] [-1;2]]
        [[0.0;0.0]  [1/3;1/3]  [-1/3;2/3]],
        [zeros(1,3); 0 1 1],
        extradimensions=["z", "sublattice"],
        highsymmetrypoints=kdict_tri
    )

    honeycomb(a::Float64=1.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  δ_hex  ],
        [ 0.0 0.0; 0 1 ];
        extradimensions=["z", "sublattice"],
        highsymmetrypoints=kdict_tri
    )

    honeycomb_AA(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  δ_hex [0.0;0.0] δ_hex],
        [ 0 0 z z; 0 1 0 1 ];
        extradimensions=["z", "sublattice"],
        highsymmetrypoints=kdict_tri
    )

    honeycomb_AB(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  δ_hex   -δ_hex   δ_hex],
        [ 0 0 z z; 0 1 0 1 ];
        extradimensions=["z", "sublattice"],
        highsymmetrypoints=kdict_tri
    )

    honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0]  -δ_hex   δ_hex   -δ_hex],
        [ 0 0 z z; 0 1 0 1 ];
        extradimensions=["z", "sublattice"],
        highsymmetrypoints=kdict_tri
    )

    # honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
    #     a .* A_tri,
    #     [[0.0;0.0]   δ_hex   δ_hex  2*δ_hex],
    #     [ 0 0 z z; 0 1 0 1 ];
    #     extradimensions=["z", "sublattice"],
    #     highsymmetrypoints=kdict_tri
    # )

    function honeycomb_twisted(N::Int, a::Float64=1.0, z::Float64=3.0)
        lat = honeycomb(a)
        return twist_triangular_2D(lat, lat, N; z=z, m=1)
    end

    function honeycomb_twisted_ABBA(N::Int, a::Float64=1.0, z::Float64=3.0)
        lat1 = honeycomb_AB(a, z)
        lat2 = honeycomb_AB(a, z)

        return twist_triangular_2D(lat1, lat2, N; z=z, m=1)
    end

    function honeycomb_twisted_ABAB(N::Int, a::Float64=1.0, z::Float64=3.0)
        lat1 = honeycomb_AB(a, z)
        lat2 = honeycomb_BA(a, z)

        return twist_triangular_2D(lat1, lat2, N; z=z, m=1)
    end


end
