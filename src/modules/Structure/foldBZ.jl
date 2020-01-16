function foldBZ!(M, kpoints::AbstractMatrix)
    """
    M = B^T B, where B=[G1 G2] contains the reciprocal lattice vectors as columns
    kpoints has the lattice points as columns in units of lattice vectors
    """
    G0sq = sum(M)
    b = M[1,2]/G0sq

    for j_ = 1:size(kpoints,2)
        kpoints[:,j_] .= mod.(kpoints[:,j_], 1.0)
        α = M * kpoints[:,j_]

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
        kpoints[:,j_] -= δk
    end

    nothing
end
foldBZ!(lat::Lattice, kpoints::AbstractMatrix) = foldBZ!(transpose(getB(lat))*getB(lat),kpoints)