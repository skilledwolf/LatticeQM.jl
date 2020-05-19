
function foldcell!(M, points::AbstractMatrix)
    @assert size(M)==(2,2) "Cell folding is only supported for d=2 lattices."

    """
    Fold all points into the two-dimensional cell with the metric M.

    M = B^T B, where B=[G1 G2] contains the reciprocal lattice vectors as columns
    points has the lattice points as columns in units of lattice vectors
    """
    G0sq = sum(M)
    b = M[1,2]/G0sq

    for j_ = 1:size(points,2)
        points[1:2,j_] .= mod.(points[1:2,j_], 1.0)
        α = M * points[1:2,j_]

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
        points[1:2,j_] -= δk
    end

    points
end


"""
Fold coordinates of k-points into the first Brillouin zone.
"""
function foldBZ!(lat::Lattice, points::AbstractMatrix)
    d = latticedim(lat)
    B = getB(lat)[1:d,1:d]
    foldcell!(transpose(B)*B, points[1:d,:])
end


"""
Fold all latticecoordinates into the first primitive unit cell.
"""
function foldPC!(lat::Lattice; shift=0.0)
    d = latticedim(lat)
    @assert d == 2 "At the moment folding is only implented for 2d lattices."

    A = getA(lat)
    lat.spacecoordinates .-= shift

    # This piece of code handles the special case of a triangular lattice.
    # We ensure that the primitive lattice vectors have angle 2π/3, not 2π/6
    # (it's equivalent, but foldcell! assumes the former)
    # Note: this part is not thoroughly tested.
    α = acos(dot(A[:,1],A[:,2])/(norm(A[:,1])*norm(A[:,2])))/(2π)
    if norm(α)≈1/6
        # print("Changing lattice basis...")
        specialpoints = deepcopy(lat.specialpoints)

        # println("Modifying lattice vectors...")
        T = [1 -1*sign(α); 0 1*sign(α)]
        lat.basis[1:2,1:2] = A[1:2,1:2] * T
        lat.spacecoordinates[1:2,:] = inv(T) * lat.spacecoordinates[1:2,:]

        for (k,v) in lat.specialpoints.coord # update high-symmetry points
            specialpoints.coord[k] = transpose(T) * v
        end
        lat.specialpoints = specialpoints
    end

    A = getA(lat)[1:d,1:d]
    foldcell!(transpose(A) * A, lat.spacecoordinates)

    lat
end

