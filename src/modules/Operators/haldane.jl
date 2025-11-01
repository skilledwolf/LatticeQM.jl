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
import ..Utils: cKDTree

function addhaldane_fast!(hops, lat::Lattice, t2::Function; ϕ=π/2, cellrange=1, mode=:none, zmode=:none)

    cross2D(x, y) = x[1] * y[2] - x[2] * y[1] # needed later on in this scope

    # Lattice references 
    neighbors = Lattices.getneighborcells(lat, cellrange; halfspace=false, innerpoints=true, excludeorigin=false)
    D=Lattices.spacedim(lat)
    N = Lattices.countorbitals(lat)
    A = Lattices.basis(lat,:, 1:Lattices.latticedim(lat))
    R0 = zero(first(neighbors))
    points = Lattices.allpositions(lat)
    poinst2 = deepcopy(points)

    # Build lookup
    trees = Dict(R => cKDTree(transpose(points[1:D,:].+A*R)) for R=neighbors)
    allpoints = hcat((points[1:D,:].+A*R for R in neighbors)...)
    largetree = cKDTree(transpose(allpoints[1:D,:]))

    # Construction
    # hops = Hops()
    for δR = neighbors
        if !haskey(hops, δR)
            hops[δR] = spzeros(ComplexF64, N, N)
        end
        if !haskey(hops, -δR)
            hops[-δR] = spzeros(ComplexF64, N, N)
        end
    end

    for (R,tree) in trees

        result = trees[R0].sparse_distance_matrix(tree, sqrt(3)+1e-3)
        filter!(x->x[2]>sqrt(3)-1e-3, result)

        # If there are no next-nearest neighbor pairs for this R, skip expensive work
        if isempty(result)
            continue
        end

        points2 = deepcopy(points)
        
        points2[1:D,:] .+= A*R
        midpoints = hcat(((points[1:D,j+1]+points2[1:D,i+1])/2 for (i,j)=keys(result))...)

        # println(R, count(map(x->length(x)==0, midids)))
        ids = map(first, largetree.query_ball_point(transpose(midpoints), 0.53))

        # hops[R] = spzeros(ComplexF64, N,N)
        for (i,(k,v)) in enumerate(result)
            r0 = allpoints[:,ids[i]+1]
            r1 = points[:,k[1]+1]
            r2 = points2[:,k[2]+1]

            x = r1[1:D]-r0[1:D]; y=r2[1:D]-r0[1:D]
            v = t2((r1+r2)/2) * exp(1.0im * ϕ * sign( cross2D(-y, x) ) )
            hops[R][Iterators.reverse(k.+1)...] += v
            if R!=R0
                hops[-R][(k.+1)...] += conj(v)
            end
        end
        # h[-R] = deepcopy(transpose(h[R]))
    end

    hops
end