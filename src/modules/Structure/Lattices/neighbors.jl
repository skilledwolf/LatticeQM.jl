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

    for (δR, δA) in zip(neighbors, δAs)
        for i=1:N
            Ri = R[:,i] .+ δA
            for j=1:N
                Rj = R[:,j]

                if d-1e-10 < norm(Ri-Rj) < d+1e-10 #&& abs2(Z[i]-Z[j]) < 0.01
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
    getneighborcells(A, k=1; halfspace=true, innerpoints=false, excludeorigin=true)

A naive implementation to find a list of `k`-th-nearest neighboring unit cells given lattice vectors `A[:,i]`.
If `halfspace=true`, the list only contain `[I,J]` without its partner `[-I,-J]`.
If `innerpoints=true`, returns all neighboring cells up to and including the `k`-th ones.
"""
function getneighborcells(A::AbstractMatrix, k::Int=1; halfspace=true, innerpoints=false, excludeorigin=true)

    ldim = size(A,2)

    if ldim == 0 # special case
        return excludeorigin ? Vector{Int64}[] : Vector{Int64}[Int64[]]
    end

    n = ceil(Int,sqrt(3)*(k+1))
    IJ  = hcat([[x...] for x = Iterators.product(Iterators.repeated(-n:n, ldim)...)]...)

    D = map(x->round(LinearAlgebra.norm(x); digits=9), eachcol(A*IJ))
    d = sort(unique(D)) # unique distances from origin unit cell

    IJ = innerpoints ? IJ[:,D .<= d[k+1]] : IJ[:,D .== d[k+1]]

    # naive iteration to make sure [I,J] and [-I,-J] do not both appear
    IJ = [Vector(v) for v=eachcol(IJ)]
    if halfspace
        if excludeorigin
            for el=IJ
                filter!(x->x != -1 .* el, IJ)
            end
        else
            for el=IJ
                filter!(x->(x==zero(el)) || (x != -1 .* el), IJ)
            end
        end
    else 
        if excludeorigin
            filter!(x->x != -1 .* x, IJ) # remove the origin
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

### THE FOLLOWING IS NOT WOKRING CORRECTLY. BEST DELETE IT...
# """
#     commonneighbor(i,j,NN)

# Given a dictionary `NN` that contains lists of nearest neighbors, we search for the
# common nearest neighbor of site index i and site index j.

# This method is primarily used in the construction of Haldane type models, where the common nearest
# neighbor determines the chirality of the next-nearest neighbor bond.
# """
# function commonneighbor(i,j,NN)

#     k = -1; R0 = zero(first(keys(NN)))
#     for (δR0, NN_pairs) = NN
#         ind1 = findall(pair->pair[1]==i, NN_pairs)
#         ind2 = findall(pair->pair[2]==j, NN_pairs)

#         if isempty(ind1) || isempty(ind2)
#             continue
#         else
#             pairs1 = NN_pairs[ind1]
#             pairs2 = NN_pairs[ind2]
#         end

#         for pair1=pairs1, pair2=pairs2
#             if pair1[2]==pair2[1]
#                 k = pair1[2]
#                 break
#             end
#         end

#         if k > -1
#             R0 = δR0
#             break
#         end
#     end

#     return k, R0
# end
