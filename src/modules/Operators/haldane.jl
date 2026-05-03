import ..Structure.Lattices

function gethaldane(args...; kwargs...)
    newhops = Hops()
    addhaldane!(newhops, args...; kwargs...)
end

"""
addhaldane!(hops, lat, t2; ϕ=π/2, spinhalf=false, cellrange=1, mode=:none, zmode=:none)

This method is a somewhat inefficient way to compute the haldane hopping matrix.
The only upside to it is that it uses methods that I already implemented and
that it is fairly general.
"""
addhaldane!(hops, lat::Lattice, t2::Number; kwargs...) = addhaldane!(hops, lat, x->t2; kwargs...)
function addhaldane!(hops, lat::Lattice, t2::Function; mode=:auto, kwargs...)
    if mode==:fast || (mode==:auto && Lattices.countorbitals(lat)>200)
        addhaldane_fast!(hops, lat, t2; kwargs...)
    else
        addhaldane_naive!(hops, lat, t2; kwargs...)
    end
end

function addhaldane_naive!(hops, lat::Lattice, t2::Function; ϕ=π/2, cellrange=1, mode=:none, zmode=:none)

    d=1    
    cross2D(x, y) = x[1] * y[2] - x[2] * y[1] # needed later on in this scope

    # NN  = find_neighbors(lat, 1.0)
    NNN = Lattices.getneighbors(lat, √3; cellrange=cellrange)

    N = Lattices.countorbitals(lat)
    D = Lattices.spacedim(lat)
    ldim = Lattices.latticedim(lat)
    R = Lattices.allpositions(lat) # positions of atoms within unit cell

    A = Lattices.getA(lat)[:,1:ldim]

    neighbors = Lattices.getneighborcells(lat, 1; halfspace=false, innerpoints=true, excludeorigin=false) #[[i;j] for i=-1:1 for j=-1:1]
    δAs = [A * v for v in neighbors]

    # hops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()
    for δR = keys(NNN)
        if !haskey(hops, δR)
            hops[δR] = spzeros(ComplexF64, N*d, N*d)
        end
    end

    for (δR, NNN_pairs) = NNN
        for (i,j) in NNN_pairs
            Ri = R[:,i]
            Ri[1:D] .+= (A * δR)
            Rj = R[:,j]

            # to do: find short cut for identifiying the common nearest neighbor, instead of iterating
            # over **all** atoms in multiple (possibly moiré) unit cells
            for δAi=δAs, δri=eachcol(R[1:D,:])
                R0 = δAi .+ δri
                
                x = Ri[1:D]-R0[1:D]; y=Rj[1:D]-R0[1:D]
                if isapprox(norm(x), 1; atol=sqrt(eps())) && isapprox(norm(y), 1; atol=sqrt(eps())) # we found the common
                    hops[δR][1+d*(i-1):d+d*(i-1), 1+d*(j-1):d+d*(j-1)] += t2((Ri+Rj)/2) * I * exp(1.0im * ϕ * sign( cross2D(-y, x) ) )#* 1im * sign( cross2D(R0.-Rj, Ri.-R0) ) #
                    break
                end
            end
        end
    end

    hops
end


import SparseArrays: spzeros
import ..TightBinding
import ..Structure
import NearestNeighbors: KDTree, inrange

# Find all pairs (i, j) such that ‖ptsA[:, i] − ptsB[:, j]‖ ≤ radius. Returns
# Dict mapping (i, j) → distance. Replaces scipy.spatial.cKDTree's
# sparse_distance_matrix without the PyCall dependency.
function pairs_within_radius(treeA::KDTree, ptsA::AbstractMatrix,
                             ptsB::AbstractMatrix, radius::Real)
    result = Dict{Tuple{Int,Int},Float64}()
    @views for j in axes(ptsB, 2)
        idxs = inrange(treeA, ptsB[:, j], radius)
        for i in idxs
            result[(i, j)] = norm(ptsA[:, i] .- ptsB[:, j])
        end
    end
    result
end

function addhaldane_fast!(hops, lat::Lattice, t2::Function; ϕ=π/2, cellrange=1, mode=:none, zmode=:none)
    cross2D(x, y) = x[1] * y[2] - x[2] * y[1]

    NNN          = Lattices.getneighbors(lat, sqrt(3); cellrange=cellrange)
    cellneighbors = Lattices.getneighborcells(lat, cellrange; halfspace=false, innerpoints=true, excludeorigin=false)
    D            = Lattices.spacedim(lat)
    N            = Lattices.countorbitals(lat)
    A            = Lattices.basis(lat, :, 1:Lattices.latticedim(lat))
    points       = Lattices.allpositions(lat)

    # Stack atom positions across all cellrange-shifted cells; one KDTree query
    # finds the common nearest neighbour of any NNN pair in O(log Ncells·N).
    # This is the only place a spatial index is asymptotically helpful — pair
    # discovery already uses Lattices.getneighbors which is itself O(N·logN).
    allpoints = hcat((points[1:D, :] .+ A * R for R in cellneighbors)...)
    largetree = KDTree(allpoints)

    for δR in keys(NNN)
        haskey(hops, δR) || (hops[δR] = spzeros(ComplexF64, N, N))
    end

    for (δR, pairs) in NNN
        AδR = A * δR
        for (i, j) in pairs
            r1 = points[1:D, i] .+ AδR
            r2 = points[1:D, j]
            mid = (r1 .+ r2) / 2

            # On a hexagonal lattice the midpoint of any NNN pair coincides
            # with their common nearest-neighbour atom. 0.53 > 0.5 to absorb
            # numerical jitter; for non-hexagonal geometries the caller
            # should fall back to addhaldane_naive!.
            cands = inrange(largetree, mid, 0.53)
            isempty(cands) && continue
            r0 = allpoints[:, first(cands)]

            x = r1 .- r0
            y = r2 .- r0
            # `t2` may read extras (e.g. `addvalley!`'s closure inspects the
            # sublattice index at R[4]). Pass the full atom position
            # including extras, matching `addhaldane_naive!`.
            mid_full = (points[:, i] .+ points[:, j]) / 2
            mid_full[1:D] .= mid
            val = t2(mid_full) * exp(1.0im * ϕ * sign(cross2D(-y, x)))
            hops[δR][i, j] += val
        end
    end

    hops
end