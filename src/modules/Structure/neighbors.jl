@legacyalias getneighbors get_neighbors
function getneighbors(lat, d=1.0)

    neighbors = [[i;j] for i=-1:1 for j=-1:1] #if i+j>=0 && i>=0]

    N = countatoms(lat)
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

                if d-0.001 < norm(Ri-Rj) < d+0.001 && abs2(Z[i]-Z[j]) < 0.1
                    # pairs[δR][i,j] = 1
                    append!(pairs[δR], [(i,j)])
                end
            end
        end
    end

    return pairs
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
