
using LinearAlgebra: I

"""
    foldcell_fromneighbors!(points, gvectors[, M])

Fold columns of `points` into the first Wigner–Seitz cell using a list of
neighbor-cell vectors `gvectors` (both in fractional coordinates). The optional
metric `M = B'B` can be supplied to define distances in a non-orthogonal basis,
where `B = [G1 G2 …]` collects the reciprocal lattice vectors as columns.

Arguments
- `points::AbstractMatrix{Float64}`: columns are coordinates to be folded
  (fractional units, i.e. in the basis of the lattice vectors).
- `gvectors::Matrix{Float64}`: columns are the neighbor-cell shift vectors
  that span the Wigner–Seitz region (fractional units).
- `M::AbstractMatrix` (optional): metric used to compute distances;
  defaults to the identity (orthonormal basis).

Returns the same `points` matrix with all columns folded in place. This is the
low‑level implementation used by higher‑level `foldcell!` methods.
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
    foldcell!(points::AbstractMatrix{Float64}, basis::Matrix{Float64})

Fold columns of `points` into the first Wigner–Seitz cell defined by the
column vectors of `basis` (direct-lattice basis), then calls
`foldcell_fromneighbors!` with the metric `basis' * basis`.

The half-plane set is built from all `{-1,0,1}` combinations of the
**Lagrange-reduced** basis (mapped back to integer coefficients of `basis`),
which contains every Voronoi-relevant vector of a 2D lattice. Using only the
single shortest neighbor shell (the old behavior) never folds along
directions longer than the shortest lattice vector — e.g. on a rectangular
`diag(1,2)` cell, points were left unfolded along the long axis.
"""
function foldcell!(points::AbstractMatrix{Float64}, basis::Matrix{Float64})

    d = size(basis, 2)
    d == 0 && return points

    _, U = _lattice_reduce(basis)
    combos = [U * [x...] for x in Iterators.product(Iterators.repeated(-1:1, d)...) if any(!=(0), x)]
    vectors = float(hcat(combos...))

    foldcell_fromneighbors!(points, vectors, transpose(basis)*basis)
end
precompile(foldcell!, (Matrix{Float64}, Matrix{Float64}))


"""
    foldBZ!(points, lat::Lattice; shift=0.0)

Fold k‑points (columns of `points`) into the first Brillouin zone of lattice
`lat`. K‑points are assumed to be in fractional (reciprocal‑basis) coordinates.
An optional `shift` vector can be subtracted before folding.
"""
foldBZ!(points, lat::Lattice; shift=0.0) = (points[:,:].-=shift; foldcell!(points, getB(lat)))


"""
    foldPC!(points, lat::Lattice; shift=0.0)

Fold lattice coordinates (columns of `points`) into the first primitive unit
cell of `lat`. Coordinates must be fractional in the direct‑lattice basis.
An optional `shift` vector can be subtracted before folding.
"""
foldPC!(points, lat::Lattice; shift=0.0) = (points[:,:].-=shift; foldcell!(points, getA(lat)))

"""
    foldPC!(lat::Lattice; shift=0.0)

In‑place variant that folds `lat.spacecoordinates` into the first primitive
unit cell of `lat`. Coordinates are interpreted as fractional. Returns the
folded coordinate submatrix view for convenience.
"""
function foldPC!(lat::Lattice; shift=0.0)
    d = latticedim(lat)

    lat.spacecoordinates[:,:] .-= shift

    points = lat.spacecoordinates[1:d,:]
    foldcell!(points, getA(lat))
    lat.spacecoordinates[1:d,:] .= points

    points
end
