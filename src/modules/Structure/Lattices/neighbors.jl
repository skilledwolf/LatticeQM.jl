function getneighbors(lat, d=1.0; cellrange::Int=1)

    # neighbors = [[i;j] for i=-1:1 for j=-1:1] #if i+j>=0 && i>=0]
    neighbors = getneighborcells(lat, cellrange; halfspace=false, innerpoints=true, excludeorigin=false)
    neighbors = map(x->[x...], neighbors)

    N = countorbitals(lat)
    ldim = latticedim(lat)
    R = positions(lat)
    # Z = extracoordinates(lat,"z")

    A = getA(lat)

    δAs = [A * v for v=neighbors]

    pairs = Dict{Vector{Int}, Vector{Tuple{Int,Int}}}()
    for δR in neighbors
        pairs[δR] = Vector{Tuple{Int,Int}}() #spzeros(Int,N,N)
    end

    # Relative window: an absolute ±1e-10 was fragile both for large
    # coordinates (moiré cells — fp error grows with position magnitude) and
    # for target distances derived from the geometry rather than typed as
    # exact literals.
    tol = max(1e-10, 1e-8 * d)
    for (δR, δA) in zip(neighbors, δAs)
        for i=1:N
            Ri = R[:,i] .+ δA
            for j=1:N
                Rj = R[:,j]

                if d-tol < norm(Ri-Rj) < d+tol #&& abs2(Z[i]-Z[j]) < 0.01
                    # pairs[δR][i,j] = 1
                    append!(pairs[δR], [(i,j)])
                end
            end
        end
    end

    return pairs
end

import LinearAlgebra

"""
    _lattice_reduce(A) → (A_red, U)

Lagrange (Gauss) reduction of the lattice basis `A` (columns are lattice
vectors): returns a reduced basis `A_red = A * U` spanning the same lattice,
with `U` unimodular (integer, |det U| = 1). For a reduced 2D basis
`‖b1‖ ≤ ‖b2‖` and `|b1·b2| ≤ ½‖b1‖²` (angle in [60°,120°]), so fixed-size
integer search boxes are sufficient to enumerate near-origin lattice points —
which is false for sheared inputs like `[1 6; 0 1]`. Only `size(A,2) == 2`
is reduced; other dimensions return `U = I` unchanged.
"""
function _lattice_reduce(A::AbstractMatrix)
    ldim = size(A, 2)
    U = Matrix{Int}(LinearAlgebra.I, ldim, ldim)
    ldim == 2 || return float(A), U

    b1 = float(A[:, 1]); b2 = float(A[:, 2])
    for _ in 1:1000   # Lagrange reduction terminates fast; bound is a safety net
        if LinearAlgebra.norm(b1) > LinearAlgebra.norm(b2)
            b1, b2 = b2, b1
            U[:, 1], U[:, 2] = U[:, 2], U[:, 1]
        end
        m = round(Int, LinearAlgebra.dot(b1, b2) / LinearAlgebra.dot(b1, b1))
        m == 0 && break
        b2 = b2 .- m .* b1
        U[:, 2] .-= m .* U[:, 1]
    end
    hcat(b1, b2), U
end

"""
    getneighborcells(A, k=1; halfspace=true, innerpoints=false, excludeorigin=true)

Find the `k`-th-nearest neighboring unit cells given lattice vectors `A[:,i]`,
returned as integer coefficient vectors in the basis `A`.
If `halfspace=true`, the list only contain `[I,J]` without its partner `[-I,-J]`.
If `innerpoints=true`, returns all neighboring cells up to and including the `k`-th ones.

The enumeration runs over a fixed box in the **Lagrange-reduced** basis (see
[`_lattice_reduce`](@ref)) and maps the results back, so skewed/non-reduced
bases (e.g. sheared supercells) find their true nearest cells; a fixed box in
the raw basis silently missed them.
"""
function getneighborcells(A::AbstractMatrix, k::Int=1; halfspace=true, innerpoints=false, excludeorigin=true)

    ldim = size(A,2)

    if ldim == 0 # special case
        return excludeorigin ? Vector{Int64}[] : Vector{Int64}[Int64[]]
    end

    A_red, U = _lattice_reduce(A)

    n = ceil(Int,sqrt(3)*(k+1))
    IJ  = hcat([[x...] for x = Iterators.product(Iterators.repeated(-n:n, ldim)...)]...)

    D = map(x->round(LinearAlgebra.norm(x); digits=9), eachcol(A_red*IJ))
    d = sort(unique(D)) # unique distances from origin unit cell

    IJ = innerpoints ? IJ[:, D .<= d[k+1]] : IJ[:, D .== d[k+1]]

    # map reduced-basis coefficients back to original-basis coefficients
    IJ = [U * Vector(v) for v=eachcol(IJ)]

    if halfspace
        # keep the first representative of each ±pair (don't mutate while iterating)
        seen = Set{Vector{Int}}()
        keep = Vector{Vector{Int}}()
        for v in IJ
            if all(iszero, v)
                excludeorigin || push!(keep, v)
            elseif !(v in seen) && !(-v in seen)
                push!(keep, v)
                push!(seen, v)
            end
        end
        IJ = keep
    else
        if excludeorigin
            filter!(x->!all(iszero, x), IJ) # remove the origin
        end
    end

    return IJ
end



"""
    getneighborcells(lat, k=1; halfspace=true, innerpoints=false, excludeorigin=true)

A naive implementation to find a list of `k`-th-nearest neighboring unit cells.
If `halfspace=true`, the list only contain `[I,J]` without its partner `[-I,-J]`.
If `innerpoints=true`, returns all neighboring cells up to and including the `k`-th ones.
"""
getneighborcells(lat, args...; kwargs...) = getneighborcells(getA(lat), args...; kwargs...)


"""
    getneighborBZ(lat, k=1; halfspace=true, innerpoints=false, excludeorigin=true)

This is the analogue of method `getneighborcells()`, except that it looks for 
nearest neighbor cells in reciprocal space.
"""
getneighborBZ(lat, args...; kwargs...) = getneighborcells(getB(lat), args...; kwargs...)

import ...Utils: padvec

function getneighbordict(lat::Lattice, cellrange::Int, d::Int)
    neighbors = getneighborcells(lat, cellrange; halfspace=true, innerpoints=true, excludeorigin=false)
    A = getA(lat)
    Dict(δL => padvec(A * δL, d) for δL in neighbors)
end

function getneighbordict(lat::Lattice, cellrange::Int)
    d = allspacedim(lat)
    # neighbors = getneighborcells(lat, cellrange; halfspace=true, innerpoints=true, excludeorigin=false)
    # A = getA(lat)
    # Dict(δL => A * δL for δL in neighbors)
    getneighbordict(lat, cellrange, d)
end