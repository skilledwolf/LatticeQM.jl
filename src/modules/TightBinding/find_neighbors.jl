
function get_neighbors(lat, d=1.0)

    neighbors = [[i;j] for i=-1:1 for j=-1:1] #if i+j>=0 && i>=0]

    N = atom_count(lat)
    R = positions(lat)
    A = get_A(lat)

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

                if 0.95*d < norm(Ri-Rj) < 1.05*d
                    # pairs[δR][i,j] = 1
                    append!(pairs[δR], [(i,j)])
                end
            end
        end
    end

    return pairs
end

function find_neighbors(lat, d=1.0)
    """
    returns neighbors at distance d in format Dict(δR => [(i1,j1), (i2,j2)...], ...)
    """

    function  get_t(d=1.0)
        function tNN(r1, r2=0.0)
            δr = r1 .- r2
            dr = norm(δr)

            if - 0.1 < dr-d < 0.1 && abs(δr[3]) < 0.1
                return 1
            else
                return 0
            end
        end
    end

    neighbors = get_hops(lat, get_t(d); format=:sparse)

    return Dict(δR => map(Tuple, findall(!iszero, mat)) for (δR, mat)=neighbors)
end

function find_common_neighbor(i,j,NN)

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
