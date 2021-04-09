# __precompile__()
"""
    Geometries2D

Provides predefined lattice objects (such as two-dimensional honeycomb lattice).

### Example
```julia
lat = Geometries2D.honeycomb_twisted(11)
plot(lat, 3; supercell=[0:1,0:1])
```
"""
module Geometries2D

    using Base, LinearAlgebra

    using ..Structure: twist
    using ..Structure.Paths
    using ..Structure.Lattices

    # Dictionary to construct Paths in k space
    kdict_sq = LabeledPoints(
        ["Γ", "X", "Y", "M"],
        [[0.0;  0.0], [1/2;  0.0], [0.0;  1/2], [1/2;  1/2]],
        ["\$\\Gamma\$", "X", "Y", "M"],
        ["Y", "Γ", "M", "X", "Γ"]
    )

    kdict_tri = LabeledPoints(
        ["Γ", "K", "M", "K'", "M2", "M3", "Γ1"],
        [[0.0;  0.0], [1/3;  2/3], [1/2;  1/2], [2/3;  1/3], [0, 1/2], [1/2,0], [1.0,-1.0]],
        ["\$\\Gamma\$", "K", "M", "K'", "M\$_2\$", "M\$_3\$", "\$\\Gamma\$"],
        ["Γ", "K", "M", "K'", "Γ", "M"]
    )

    kdict_tri_mini = LabeledPoints(
        ["γ", "κ", "μ", "κ'", "μ2", "μ3", "γ1"],
        [[0.0;  0.0], [1/3;  2/3], [1/2;  1/2], [2/3;  1/3], [0, 1/2], [1/2,0], [1.0,-1.0]],
        ["\$\\gamma\$", "\$\\kappa\$", "\$\\mu\$", "\$\\kappa'\$", "\$\\mu_2\$", "\$\\mu_3\$", "\$\\gamma\$"],
        ["γ", "κ", "μ", "κ'", "γ", "μ"]
    )

    # Should implement this (or a similar) form in the future:
    # kdict_tri = LabeledPoints(
    #     "γ" => (label="\$\\gamma\$", coord=[0.0;  0.0]),
    #     "κ" => (label="\$\\kappa\$", coord=[1/3;  2/3]),
    #     "μ" => (label="\$\\mu\$", coord=[1/2;  1/2]),
    #     "κ'" => (label="\$\\kappa'\$", coord=[2/3;  1/3]),
    # )

    square(a::Float64=1.0) = Lattice(
        [[a,0,0] [0,a,0]],
        zeros(3,1);
        extralabels=String[],
        specialpoints=kdict_sq
    )

    #### Define simple honeycomb systems
    A_tri = [[cos(pi/6);  -sin(pi/6); 0]  [cos(pi/6);  sin(pi/6); 0]]
    A_hex = sqrt(3) .* A_tri
    δ_hex = [1/3; 1/3; 0]

    triangular(a::Float64=1.0) = Lattice(
        a .* A_tri,
        zeros(3,1);
        specialpoints=kdict_tri
    )

    # [[cos(2*pi/6);  -sin(2*pi/6)]  [cos(2*pi/6);  sin(2*pi/6)]]
    triangular_supercell(a::Float64=1.0) = Lattice(
        a .* A_tri * [[1;1] [-1;2]], #[[2;-1] [-1;2]]
        [[0.0;0.0;0.0]  [1/3;1/3;0.0]  [-1/3;2/3;0.0]],
        1.0*hcat([1 2 3]),
        extralabels=["sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb(a::Float64=1.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0;0.0]  δ_hex  ],
        1.0*hcat([0 1]);
        extralabels=["sublattice"],
        specialpoints=kdict_tri
    )

    graphene() = honeycomb(1.42)

    honeycomb_bilayer(a::Float64=1.0, z::Float64=3.0; δ::Vector{Float64}=[0.0,0.0]) = Lattice(
        a .* A_hex,
        [[0.0;0.0;-z/2]  [1/3; 1/3; -z/2] [δ; z/2] δ_hex+[δ; z/2]],
        1.0*hcat([0 1 0 1]);
        extralabels=["sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb_AA(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0;-z/2]  [1/3; 1/3; -z/2] [0.0;0.0;z/2] [1/3; 1/3; z/2]],
        1.0*hcat([0 1 0 1]);
        extralabels=["sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb_AB(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0;-z/2]  [1/3; 1/3; -z/2]  [2/3; 2/3; z/2]   [1/3; 1/3; z/2]],
        1.0*hcat([0 1 0 1]);
        extralabels=["sublattice"],
        specialpoints=kdict_tri
    )

    honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
        a .* A_hex,
        [[0.0;0.0;-z/2]  [-1/3; -1/3; -z/2]   [1/3; 1/3; z/2]   [-1/3; -1/3; z/2]],
        1.0*hcat([0 1 0 1]);
        extralabels=["sublattice"],
        specialpoints=kdict_tri
    )

    # honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
    #     a .* A_tri,
    #     [[0.0;0.0]   δ_hex   δ_hex  2*δ_hex],
    #     [ 0 0 z z; 0 1 0 1 ];
    #     extralabels=["z", "sublattice"],
    #     specialpoints=kdict_tri
    # )

    import ..Structure

    function honeycomb_twisted(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
        lat = honeycomb(a)
        slat = twist(lat, lat, N; z=z, m=1)
        slat.specialpoints = kdict_tri_mini
        if fold 
            Structure.Lattices.foldPC!(slat; shift=[1/3,1/3,0])
        end
        return slat
    end

    function honeycomb_twisted_ABBA(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
        lat1 = honeycomb_AB(a, z)
        translate!(lat1, 3, z/2)
        lat2 = honeycomb_AB(a, z)
        translate!(lat2, 3, z/2)
        slat = twist(lat1, lat2, N; z=z, m=1)
        slat.specialpoints = kdict_tri_mini
        if fold 
            Structure.Lattices.foldPC!(slat; shift=[1/3,1/3,0])
        end
        return slat
    end

    function honeycomb_twisted_ABAB(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
        lat1 = honeycomb_AB(a, z)
        translate!(lat1, 3, z/2)
        lat2 = honeycomb_BA(a, z)
        translate!(lat2, 3, z/2)
        slat = twist(lat1, lat2, N; z=z, m=1)
        slat.specialpoints = kdict_tri_mini
        if fold 
            Structure.Lattices.foldPC!(slat; shift=[1/3,1/3,0])
        end
        return slat
    end


    using ..Structure.Lattices: displaceZ!, getneighborBZ

    function smoothdisplaceZ!(lat, δz_even=0.055, δz_odd=0.0; sharp::Real=1)
        @assert latticedim(lat) == 2 "Lattice must be two-dimensional."
        @assert abs(dot(getA(lat,:,1), getA(lat,:,2))/(norm(getA(lat,:,1))*norm(getA(lat,:,2)))) ≈ 0.5

        # Define smooth (super)lattice functions
        Gs = [getB(lat)*v for v = getneighborBZ(lat,1; halfspace=false, innerpoints=true)]

        anglef(b1,b2) = acos(dot(b1,b2)/norm(b1)/norm(b2))
        angles = [anglef(Gs[1], G) for G=Gs]
        angles /= minimum(abs.(filter(α->α!=0, angles)))
        signs = (-1).^round.(Int,angles)
        
        fodd(x::AbstractVector) =  1/(3*√3) * sum((s/2)*sin(2*π*dot(G,x)) for (s,G)=zip(signs,Gs))
        foddsharp(x::AbstractVector) = tanh(2 * sharp * fodd(x))
        feven(x::AbstractVector) = 1/3 + 2/9 * sum((1/2)*cos(2*π*dot(G,x)) for G=Gs)
        fevensharp(x::AbstractVector) = tanh(sharp * feven(x))/2

        if sharp > 1
            displaceZ!(lat, p -> sign(p[3]) * (δz_even * (fevensharp(p)-0.5) + δz_odd * foddsharp(p)))
        else
            displaceZ!(lat, p -> sign(p[3]) * (δz_even * (feven(p)-0.5) + δz_odd * fodd(p)))
        end
    end

end
