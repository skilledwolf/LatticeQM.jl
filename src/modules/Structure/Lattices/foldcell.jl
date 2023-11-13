
using LinearAlgebra: I

"""
    foldcell_fromneighbors!(points::Matrix{Float64}, vectors::Matrix{Float64}, M=id)

Fold all points into the first Wigner-Seitz cell, given discrete points `points`, lattice vectors `vectors`
(both in fracdtional coordinates) pointing to all neighboring unit cells, and given the lattice metric M=B^T*B.

M = B^T B, where B=[G1 G2] contains the (reciprocal) lattice vectors as columns
points has the lattice points as columns in units of lattice vectors.

This should be considered the low-level API for foldcell! methods, and implements a general folding algorithm.
"""
foldcell_fromneighbors!(points::AbstractMatrix{Float64}, gvectors::Matrix{Float64}) = foldcell_fromneighbors!(points, gvectors, 1.0*I)
function foldcell_fromneighbors!(points::AbstractMatrix{Float64}, gvectors::Matrix{Float64}, M::AbstractMatrix)


    norms = [dot(v,M,v) for v=eachcol(gvectors)]
    distance(i, p) = (g=gvectors[:,i]; N=norms[i]; dot((p-g/2), M, g)/N)

    shift = true

    while shift

        shift = false

        for (i, g) in enumerate(eachcol(gvectors))
            for (j,p) in enumerate(eachcol(points))
                d = distance(i,p)

                if d > 0
                    points[:,j] .-= ceil(d) * g
                    shift = true # we had to shift an atom. we need to do one more loop to find out if it needs more shifts
                end
            end
        end

    end 

    points
end
precompile(foldcell_fromneighbors!, (Matrix{Float64}, Matrix{Float64}, Matrix{Float64}))



"""
    foldcell!(points::Matrix{Float64}, basis::Matrix{Float64})

Same `foldcell_fromneighbors!` but constructs the nearest-neighboring cells for a given lattice basis `basis`
using `getneighborcells`.

"""
function foldcell!(points::AbstractMatrix{Float64}, basis::Matrix{Float64})

    vectors = getneighborcells(basis, 1; halfspace=false, innerpoints=true, excludeorigin=true)
    vectors = float(hcat(vectors...))

    foldcell_fromneighbors!(points, vectors, transpose(basis)*basis)
end
precompile(foldcell!, (Matrix{Float64}, Matrix{Float64}))


"""
    foldBZ!(points, lat::Lattice; shift=0.0)

Fold coordinates of k-points into the first Brillouin zone. Note that k-points are assumed to
be in fractional coordinates.
"""
foldBZ!(points, lat::Lattice; shift=0.0) = (points[:,:].-=shift; foldcell!(points, getB(lat)))


"""
    foldPC!(points, lat::Lattice; shift=0.0)

Fold all fractional lattice coordinates `points` into the first primitive unit cell.
"""
foldPC!(points, lat::Lattice; shift=0.0) = (points[:,:].-=shift; foldcell!(points, getA(lat)))

"""
    foldPC!(lat::Lattice; shift=0.0)

Fold all fractional lattice coordinates `lat.spacecoordinates` into the first primitive unit cell.
"""
function foldPC!(lat::Lattice; shift=0.0)
    d = latticedim(lat)

    lat.spacecoordinates[:,:] .-= shift

    points = lat.spacecoordinates[1:d,:]
    foldcell!(points, getA(lat))
    lat.spacecoordinates[1:d,:] .= points

    points
end

