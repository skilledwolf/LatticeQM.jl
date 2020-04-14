@legacyalias getneighbors get_neighbors
function getneighbors(lat, d=1.0; cellrange::Int=1)

    # neighbors = [[i;j] for i=-1:1 for j=-1:1] #if i+j>=0 && i>=0]
    neighbors = getneighborcells(lat, cellrange; halfspace=false, innerpoints=true, excludeorigin=false)
    neighbors = map(x->[x...], neighbors)

    N = countorbitals(lat)
    R = positions(lat)
    Z = extrapositions(lat,"z")

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

                if d-0.01 < norm(Ri-Rj) < d+0.01 && abs2(Z[i]-Z[j]) < 0.01
                    # pairs[δR][i,j] = 1
                    append!(pairs[δR], [(i,j)])
                end
            end
        end
    end

    return pairs
end


"""
    A naive implementation to find a list of n-th-nearest neighboring unit cells.
    If halfspace=true, the list only contain [I,J] without its partner [-I,-J].
"""
function getneighborcells(lat, k::Int=1; halfspace=true, innerpoints=false, excludeorigin=true)


    ldim = latticedim(lat)
    A = getA(lat)

    if ldim == 0 # special case
        return excludeorigin ? Vector{Int64}[] : Vector{Int64}[Int64[]]
    end

    n = ceil(Int,sqrt(3)*(k+1))
    IJ  = collect(Iterators.product(Iterators.repeated(-n:n, ldim)...))
    IJ = map(x->[x...], IJ)

    D = map(x->round.(norm(A*[x...]); digits=7), IJ) # cell distances

    d = sort(unique(D)) # unique distances from origin unit cell

    if innerpoints
        result = IJ[D .<= d[k+1]]
    else
        result = IJ[D .== d[k+1]]
    end

    # naive iteration to make sure [I,J] and [-I,-J] do not both appear
    if halfspace
        if excludeorigin
            for el=result
                filter!(x->x != -1 .* el, result)
            end
        else
            for el=result
                filter!(x->(x==zeros(eltype(el), length(el))) || (x != -1 .* el), result)
            end
        end

    end

    return result
end

@legacyalias commonneighbor find_common_neighbor
function commonneighbor(i,j,NN)

    k = -1; R0 = zero(first(keys(NN)))
    for (δR0, NN_pairs) = NN
        ind1 = findall(pair->pair[1]==i, NN_pairs)
        ind2 = findall(pair->pair[2]==j, NN_pairs)

        if isempty(ind1) || isempty(ind2)
            continue
        else
            pairs1 = NN_pairs[ind1]
            pairs2 = NN_pairs[ind2]
        end

        for pair1=pairs1, pair2=pairs2
            if pair1[2]==pair2[1]
                k = pair1[2]
                break
            end
        end

        if k > -1
            R0 = δR0
            break
        end
    end

    return k, R0

end
