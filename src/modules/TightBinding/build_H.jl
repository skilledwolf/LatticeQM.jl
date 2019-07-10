function build_H(lat::Lattice, t::Function; mode=:nospin, format=:auto)

    build_H_sparse(get_A_3D(lat), positions3D(lat), t, mode=mode, format=format)
end

@fastmath function build_H_sparse(A::Matrix{Float64}, R::Matrix{Float64}, t::Function; mode=:nospin, format=:auto, precision::Float64=1e-8)# where {T<:AbstractMatrix{Float64}}

    # Get lattice neighbors
    neighbors = [[i;j] for i=-1:1 for j=-1:1 if i+j>=0 && i>=0]
    δA = [A*v for v in neighbors]

    hops = build_hops_sparse(R, neighbors, δA, t; precision=precision)

    # Create the Hermitian conjugates
    for (δL, T) in hops
        hops[-δL] = T'
    end

    # Build the Bloch matrix by adding the hopping matrices with the correct phases
    build_BlochH(hops; mode=mode, format=format)
end

function build_hops_sparse(R, neighbors, δA, t; precision=precision)

    N = size(R)[2]

    hops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()

    # Preallocate memory
    # Arbitrary limit for dense allocation
    if N > 500
        maxind = round(Int, 0.60 * N^2)
    else
        maxind = N^2
    end

    IS = ones(Int, maxind)
    JS = ones(Int, maxind)
    VS = zeros(ComplexF64, maxind)
    R0 = zero(R)
    V =  zeros(ComplexF64, N)

    for (δL,a) in zip(neighbors,δA)
        count = hopmat_from_gen!(IS, JS, VS, R0, V, N, maxind, R, a, t; precision=precision)
        if count > maxind
            error("Not enough memory-reserved. Your problem does not appear to be sparse.")
        end

        hops[δL] = sparse(IS[1:count], JS[1:count], VS[1:count], N, N)
    end

    hops
end

function hopmat_from_gen!(IS::Vector{Int}, JS::Vector{Int}, VS::Vector{ComplexF64}, R0::Matrix{Float64}, V::Vector{ComplexF64}, N::Int, maxind::Int, R::Matrix{Float64}, a::Vector{Float64}, t::F; precision::Float64) where {F<:Function}

    count = 1

    R0 .= R.+a

    for j=1:N
        V .= t.(eachcol(R0.-R[:,j]))

        for i=1:N
            if abs(V[i]) > precision
                if count > maxind
                    break
                end
                IS[count] = i
                JS[count] = j
                VS[count] = V[i]

                count = count+1
            end
        end
    end

    return count-1
end
