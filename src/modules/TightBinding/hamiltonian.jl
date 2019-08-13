@fastmath function build_H(lat::Lattice, t::Function; mode=:nospin, format=:auto, precision::Float64=1e-8)# where {T<:AbstractMatrix{Float64}}

    # Get lattice neighbors
    neighbors = [[i;j] for i=-1:1 for j=-1:1 if i+j>=0 && i>=0]

    hops = get_hops(lat, neighbors, t; precision=precision, format=format)

    # Build the Bloch matrix by adding the hopping matrices with the correct phases
    build_BlochH(hops; mode=mode)
end

function decide_type(hops::Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}, format)

    if format==:auto
        N = size(first(values(hops)), 1)
        if N < 301
            format=:dense
        end
    end

    if format==:dense
        hops = Dict(δL => Matrix(t) for (δL, t) in hops)
    end

    hops
end

function get_hops(lat::Lattice, neighbors::Vector{Vector{Int}}, t::Function; kwargs...)
    get_hops(get_A_3D(lat), positions3D(lat), neighbors, t; kwargs...)
end

function get_hops(A::Matrix{Float64}, R::Matrix{Float64}, neighbors::Vector{Vector{Int}}, t::Function; precision::Float64=1e-8, format=:auto)

    # Get lattice neighbors
    δA = [A*v for v in neighbors]

    hops = build_hopmats(R, neighbors, δA, t; precision=precision)

    # Create the Hermitian conjugates
    for (δL, T) in hops
        hops[-δL] = T'
    end

    hops = decide_type(hops, format)

    hops
end

function build_hopmats(R, neighbors, δA, t; precision=precision)

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

    # Preallocate memory: important for huge sparse matrices
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
