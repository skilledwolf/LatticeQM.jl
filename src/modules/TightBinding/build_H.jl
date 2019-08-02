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
    # Semi-arbitrary limit for dense allocation
    if N > 500
        maxind = round(Int, 0.60 * N^2)
    else
        maxind = N^2
    end

    v0 = t(R[:,1])
    if isa(v0, Number)
        d0=1
    elseif isa(v0, AbstractMatrix)
        d0=size(v0)[1]
    else
        error("Hopping function t must return number or matrix of number")
    end

    IS = ones(Int, maxind)
    JS = ones(Int, maxind)
    VS = zeros(ComplexF64, maxind*d0, d0)
    R0 = zero(R)
    V =  zeros(ComplexF64, N*d0, d0)

    for (δL,a) in zip(neighbors,δA)
        count = hopmat_from_gen!(IS, JS, VS, R0, V, N, maxind, R, a, t, d0; precision=precision)
        if count > maxind
            error("Not enough memory-reserved. Your problem does not appear to be sparse enough.")
        end

        hops[δL] = sparse(IS[1:count], JS[1:count], VS[1:count], N, N)
    end

    hops
end

function hopmat_from_gen!(IS::Vector{Int}, JS::Vector{Int}, VS::Array{ComplexF64}, R0::Matrix{Float64}, V::Array{ComplexF64}, N::Int, maxind::Int, R::Matrix{Float64}, a::Vector{Float64}, t::F, d::Int; precision::Float64) where {F<:Function}

    count = 1

    R0 .= R.+a

    for j=1:N
        V .= vcat(t.(eachcol(R0.-R[:,j]))...)

        for i=1:N
            if abs(V[i]) > precision
                if count > maxind
                    break
                end
                IS[count] = i
                JS[count] = j
                VS[1+(count-1)*d:count*d,1:d] .= V[1+(i-1)*d:i*d, 1:d]

                count = count+1
            end
        end
    end

    return count-1
end
