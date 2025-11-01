
import ..Paths
import ..Lattices
import ..Lattices: Lattice

# Dictionary to construct Paths in k space
kdict_sq = Paths.LabeledPoints(
    ["Γ", "X", "Y", "M"],
    [[0.0;  0.0], [1/2;  0.0], [0.0;  1/2], [1/2;  1/2]],
    ["\$\\Gamma\$", "X", "Y", "M"],
    ["Y", "Γ", "M", "X", "Γ"]
)

kdict_tri = Paths.LabeledPoints(
    ["Γ", "K", "M", "K'", "M2", "M3", "Γ1"],
    [[0.0;  0.0], [1/3;  2/3], [1/2;  1/2], [2/3;  1/3], [0, 1/2], [1/2,0], [1.0,-1.0]],
    ["\$\\Gamma\$", "K", "M", "K'", "M\$_2\$", "M\$_3\$", "\$\\Gamma\$"],
    ["Γ", "K", "M", "K'", "Γ", "M"]
)

kdict_tri_mini = Paths.LabeledPoints(
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

"""
    square(a=1.0)

Square Bravais lattice with lattice constant `a`. Returns a `Lattice` with one
orbital per cell and standard special k‑points (Γ, X, Y, M).
"""
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

"""
    triangular(a=1.0)

Triangular Bravais lattice with lattice constant `a`. Returns a `Lattice` with
one orbital per cell and triangular‑lattice special points (Γ, K, M, …).
"""
triangular(a::Float64=1.0) = Lattice(
    a .* A_tri,
    zeros(3,1);
    specialpoints=kdict_tri
)

# [[cos(2*pi/6);  -sin(2*pi/6)]  [cos(2*pi/6);  sin(2*pi/6)]]
"""
    triangular_supercell(a=1.0)

A 3‑site triangular supercell useful for toy models (e.g. Kekulé). Provides a
`"sublattice"` extra label with values 1,2,3.
"""
triangular_supercell(a::Float64=1.0) = Lattice(
    a .* A_tri * [[1;1] [-1;2]], #[[2;-1] [-1;2]]
    [[0.0;0.0;0.0]  [1/3;1/3;0.0]  [-1/3;2/3;0.0]],
    1.0*hcat([1 2 3]),
    extralabels=["sublattice"],
    specialpoints=kdict_tri
)

"""
    honeycomb(a=1.0)

Two‑site honeycomb (graphene) lattice with lattice constant `a` (distance
between nearest neighbors equals `a`). Adds a `"sublattice"` extra coordinate
with values 0 (A) and 1 (B).
"""
honeycomb(a::Float64=1.0) = Lattice(
    a .* A_hex,
    [[0.0;0.0;0.0]  δ_hex  ],
    1.0*hcat([0 1]);
    extralabels=["sublattice"],
    specialpoints=kdict_tri
)

"""
    graphene()

Convenience alias for `honeycomb(1.42)` using a typical C–C bond length
in Ångström units.
"""
graphene() = honeycomb(1.42)

"""
    honeycomb_bilayer(a=1.0, z=3.0; δ=[0.0, 0.0])

Bernal‑stacked honeycomb bilayer (AB) with interlayer distance `z`. Optional
in‑plane shift `δ` (fractional coordinates) moves the top layer prior to
stacking. Extra coordinate `"sublattice"` is provided.
"""
honeycomb_bilayer(a::Float64=1.0, z::Float64=3.0; δ::Vector{Float64}=[0.0,0.0]) = Lattice(
    a .* A_hex,
    [[0.0;0.0;-z/2]  [1/3; 1/3; -z/2] [δ; z/2] δ_hex+[δ; z/2]],
    1.0*hcat([0 1 0 1]);
    extralabels=["sublattice"],
    specialpoints=kdict_tri
)

"""
    honeycomb_AA(a=1.0, z=3.0)

AA‑stacked honeycomb bilayer with interlayer distance `z`.
"""
honeycomb_AA(a::Float64=1.0, z::Float64=3.0) = Lattice(
    a .* A_hex,
    [[0.0;0.0;-z/2]  [1/3; 1/3; -z/2] [0.0;0.0;z/2] [1/3; 1/3; z/2]],
    1.0*hcat([0 1 0 1]);
    extralabels=["sublattice"],
    specialpoints=kdict_tri
)

"""
    honeycomb_AB(a=1.0, z=3.0)

AB (Bernal) stacking for a honeycomb bilayer with interlayer distance `z`.
"""
honeycomb_AB(a::Float64=1.0, z::Float64=3.0) = Lattice(
    a .* A_hex,
    [[0.0;0.0;-z/2]  [1/3; 1/3; -z/2]  [2/3; 2/3; z/2]   [1/3; 1/3; z/2]],
    1.0*hcat([0 1 1 0]);
    extralabels=["sublattice"],
    specialpoints=kdict_tri
)

"""
    honeycomb_BA(a=1.0, z=3.0)

BA stacking variant for a honeycomb bilayer with interlayer distance `z`.
"""
honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
    a .* A_hex,
    [[0.0;0.0;-z/2]  [-1/3; -1/3; -z/2]   [1/3; 1/3; z/2]   [-1/3; -1/3; z/2]],
    1.0*hcat([0 1 0 1]);
    extralabels=["sublattice"],
    specialpoints=kdict_tri
)

"""
    honeycomb_ABC(a=1.0, z=3.0)

ABC (rhombohedral) stacked tri‑layer honeycomb. Adds `"sublattice"` and
`"layer"` extra coordinates.
"""
honeycomb_ABC(a::Float64=1.0, z::Float64=3.0) = Lattice(
    a .* A_hex,
    [[0.0;0.0;0.0]  [1/3; 1/3; 0.0]  [1/3; 1/3; z]   [2/3; 2/3; z]   [-1/3; -1/3; 2*z]   [0; 0; 2*z]],
    1.0*[[0 1 0 1 0 1]; [-1 -1 1 1 -1 -1]];
    extralabels=["sublattice", "layer"],
    specialpoints=kdict_tri
)

# honeycomb_BA(a::Float64=1.0, z::Float64=3.0) = Lattice(
#     a .* A_tri,
#     [[0.0;0.0]   δ_hex   δ_hex  2*δ_hex],
#     [ 0 0 z z; 0 1 0 1 ];
#     extralabels=["z", "sublattice"],
#     specialpoints=kdict_tri
# )

"""
    triangular_twisted(N, a=1.0, z=3.0; fold=true)

Commensurate twisted bilayer built from triangular lattices using twist index
`N` (moiré periodicity). If `fold=true`, folds positions back to the primitive
cell. Sets a compact set of special k‑points.
"""
function triangular_twisted(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
    lat = triangular(a)
    slat = Lattices.twist(lat, lat, N; z=z, m=1)
    slat.specialpoints = kdict_tri_mini
    if fold 
        Lattices.foldPC!(slat; shift=[1/3,1/3,0])
    end
    return slat
end
precompile(triangular_twisted, (Int, Float64, Float64))

"""
    honeycomb_twisted(N, a=1.0, z=3.0; fold=true)

Commensurate twisted bilayer graphene (TBG) with twist index `N`. If `fold=true`
the structure is wrapped to the primitive cell.
"""
function honeycomb_twisted(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
    lat = honeycomb(a)
    slat = Lattices.twist(lat, lat, N; z=z, m=1)
    slat.specialpoints = kdict_tri_mini
    if fold 
        Lattices.foldPC!(slat; shift=[1/3,1/3,0])
    end
    return slat
end
precompile(honeycomb_twisted, (Int, Float64, Float64))

"""
    honeycomb_twisted_ABBA(N, a=1.0, z=3.0; fold=true)

Four‑layer twisted stack with ABBA stacking within layers prior to twisting.
"""
function honeycomb_twisted_ABBA(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
    lat1 = honeycomb_AB(a, z)
    Lattices.translate!(lat1, 3, z/2)
    lat2 = honeycomb_AB(a, z)
    Lattices.translate!(lat2, 3, z/2)
    slat = Lattices.twist(lat1, lat2, N; z=z, m=1)
    slat.specialpoints = kdict_tri_mini
    if fold 
        Lattices.foldPC!(slat; shift=[1/3,1/3,0])
    end
    return slat
end
precompile(honeycomb_twisted_ABBA, (Int, Float64, Float64))

"""
    honeycomb_twisted_ABAB(N, a=1.0, z=3.0; fold=true)

Four‑layer twisted stack with ABAB stacking within layers prior to twisting.
"""
function honeycomb_twisted_ABAB(N::Int, a::Float64=1.0, z::Float64=3.0; fold=true)
    lat1 = honeycomb_AB(a, z)
    Lattices.translate!(lat1, 3, z/2)
    lat2 = honeycomb_BA(a, z)
    Lattices.translate!(lat2, 3, z/2)
    slat = Lattices.twist(lat1, lat2, N; z=z, m=1)
    slat.specialpoints = kdict_tri_mini
    if fold 
        Lattices.foldPC!(slat; shift=[1/3,1/3,0])
    end
    return slat
end
precompile(honeycomb_twisted_ABAB, (Int, Float64, Float64))


import LinearAlgebra: norm, dot
import AngleBetweenVectors

"""
    smoothdisplaceZ!(lat, δz_even=0.055, δz_odd=0.0; sharp=1)

Apply a smooth out‑of‑plane displacement pattern to a hexagonal bilayer `lat`
parameterized by even/odd amplitudes. The optional `sharp` parameter controls
the smoothness (higher values approach a sign‑like modulation).
"""
function smoothdisplaceZ!(lat, δz_even=0.055, δz_odd=0.0; sharp::Real=1)
    @assert Lattices.latticedim(lat) == 2 "Lattice must be two-dimensional."
    @assert (a = AngleBetweenVectors.angle(Lattices.getA(lat,:,1), Lattices.getA(lat,:,2)); a≈pi/3 || a≈2pi/3)

    # Define smooth (super)lattice functions
    Gs = [Lattices.getB(lat)*v for v = Lattices.getneighborBZ(lat,1; halfspace=false, innerpoints=true)]

    # anglef(b1,b2) = acos(dot(b1/norm(b1),b2/norm(b2)))
    angles = [AngleBetweenVectors.angle(Gs[1], G) for G=Gs] 
    angles /= minimum(abs.(filter(α->α!=0, angles)))
    signs = (-1).^round.(Int,angles)
    
    fodd(x::AbstractVector) =  1/(3*√3) * sum((s/2)*sin(2*π*dot(G,x)) for (s,G)=zip(signs,Gs))
    foddsharp(x::AbstractVector) = tanh(2 * sharp * fodd(x))
    feven(x::AbstractVector) = 1/3 + 2/9 * sum((1/2)*cos(2*π*dot(G,x)) for G=Gs)
    fevensharp(x::AbstractVector) = tanh(sharp * feven(x))/2

    if sharp > 1
        Lattices.displaceZ!(lat, p -> sign(p[3]) * (δz_even * (fevensharp(p)-0.5) + δz_odd * foddsharp(p)))
    else
        Lattices.displaceZ!(lat, p -> sign(p[3]) * (δz_even * (feven(p)-0.5) + δz_odd * fodd(p)))
    end
end
