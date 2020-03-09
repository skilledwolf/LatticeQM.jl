
@legacyalias foldcell! foldBZ!
function foldcell!(M, points::AbstractMatrix)
    @assert latticedim(lat) == 2 # only implemented for 2d lattices at the moment

    """
    Fold all points into the two-dimensional cell with the metric M.

    M = B^T B, where B=[G1 G2] contains the reciprocal lattice vectors as columns
    points has the lattice points as columns in units of lattice vectors
    """
    G0sq = sum(M)
    b = M[1,2]/G0sq

    for j_ = 1:size(points,2)
        points[:,j_] .= mod.(points[:,j_], 1.0)
        α = M * points[:,j_]

        α1 = α[1]/M[1,1]
        α2 = α[2]/M[2,2]
        α3 = (α[1]+α[2])/G0sq

        if α1 > 1/2 && α2 < 1/2+b
            δk = [1.0, 0.0]
        elseif α1 < 1/2+b && α2 > 1/2
            δk = [0.0, 1.0]
        elseif α3 > 1/2
            δk = [1.0, 1.0]
        else
            δk = [0.0,0.0]
        end
        points[:,j_] -= δk
    end

    nothing
end

"""
Fold coordinates of k-points into the first Brillouin zone.
"""
foldBZ!(lat::Lattice, points::AbstractMatrix) = foldBZ!(transpose(getB(lat))*getB(lat),points)

"""
Fold all orbitalcoordinates into the first primitive unit cell.
"""
function foldPZ!(lat::Lattice)
    @assert latticedim(lat) == 2

    A = getA(lat)

    # This piece of code handles the special case of a triangular lattice.
    # We ensure that the primitive lattice vectors have angle 2π/3, not 2π/6
    # (it's equivalent, but foldcell! assumes the former)
    # Note: this part is not thorough tested.
    α = acos(dot(A[:,1],A[:,2])/(norm(A[:,1])*norm(A[:,2])))/(2π)
    if norm(α)≈1/6
        T = [1 -1*sign(α); 0 1*\sign(α)]
        lat.latticevectors = lat.latticevectors * T
        lat.orbitalcoordinates = inv(T) * lat.orbitalcoordinates
    end

    M = transpose(A) * A
    foldcell!(transpose(A) * A, lat.orbitalcoordinates)

    lat
end

