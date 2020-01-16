gethamiltonian(args...; mode=:nospin, kwargs...) = getbloch(gethops(args...; kwargs...); mode=mode)

###############################################################################
# Wrapper for custom types to gethops(...)
###############################################################################

addhops!(hops::AnyHops, lat::Lattice, t::Function; kwargs...) = addhops!(hops, gethops(lat, t; kwargs...))

function gethops(lat::Lattice, t::Function; format=:auto, precision::Float64=1e-8)# where {T<:AbstractMatrix{Float64}}

    # Get lattice neighbors
    neighbors = [[i;j] for i=-1:1 for j=-1:1 if i+j>=0 && i>=0]

    gethops(lat, neighbors, t; precision=precision, format=format)
end


function gethops(lat::Lattice, neighbors::Vector{Vector{Int}}, t::Function; kwargs...)
    R = allpositions(lat)
    d = size(R,1)

    A = getA(lat)
    neighbor_dict = Dict(δL => padvec(A*δL,d) for δL in neighbors)

    gethops(R, neighbor_dict, t; kwargs...)
end


function padvec(v::T, d::Int) where T<:AbstractVector{<:AbstractFloat}
"""
Make sure Vector v has length d, pad with zeros if needed.
"""
    L = length(v)
    if L > d
        error("Vector exceeds specified length $d.")
    elseif L==d
        return v
    end
    return vcat(v,zeros(d-L))
end

###############################################################################
# Main routines for gethops(...)
###############################################################################

function gethops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; precision::Float64=1e-8, format=:auto)

    N = size(R,2)
    hops = Hops() # Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()
    maxind = (N>500) ? round(Int, 0.60 * N^2) : N^2 # Preallocate memory: Semi-arbitrary limit for dense allocation

    v0 = t(R[:,1])
    if isa(v0, Number)
        d0=1
    elseif isa(v0, AbstractMatrix)
        d0=size(v0,1)
    else
        error("Hopping function t must return number or matrix of number")
    end

    # Preallocate memory: important for huge sparse matrices
    IS = Vector{Int}(undef, maxind*d0^2)
    JS = Vector{Int}(undef, maxind*d0^2)
    VS = Matrix{ComplexF64}(undef, (maxind*d0^2, 1))
    R0 = Matrix{Float64}(undef, size(R))
    V  = Matrix{ComplexF64}(undef, (N*d0, d0))

    # Heavy lifting: hopping matrix construction
    for (δL,δa) in neighbors
        R0 .= R.+δa
        hops[δL] = hoppingmatrix!(IS, JS, VS, V, R, R0, t; precision=precision)
    end

    # Create the Hermitian conjugates
    for (δL, T) in hops
        hops[-δL] = T'
    end

    decidetype(hops, format)
end


function hoppingmatrix!(IS::Vector{Int}, JS::Vector{Int}, VS::Array{ComplexF64}, V::Array{ComplexF64}, R::Matrix{Float64}, R0::Matrix{Float64}, t::F; precision::Float64) where {F<:Function}

    # Infer dimensions from array sizes
    d = size(V, 2)
    Nd = size(V,1)
    N = div(Nd,d)
    maxind = div(length(IS),d^2)

    count = 1

    for j=1:N
        tj(x) = t(x, R[:,j])

        for k=1:N # Evaluate and save the hopping elements
            V[1+d*(k-1):d+d*(k-1), 1:d] .= tj(R0[:,k])
        end

        for i=1:N
            for i0=1:d, j0=1:d
                if abs(V[(i-1)*d + i0, j0]) > precision
                    if count > maxind
                        break
                    end
                    IS[count] = (i-1)*d + i0
                    JS[count] = (j-1)*d + j0

                    VS[count] = V[(i-1)*d + i0, j0]

                    count = count+1
                end
            end
        end
    end
    count = count - 1

    if count > maxind
        error("Not enough memory-reserved. Your problem does not appear to be sparse enough.")
    end

    sparse(IS[1:count], JS[1:count], VS[1:count], N*d, N*d)
end


###################################################################################################
# Backwards compatibility
###################################################################################################
@legacyalias gethops get_hops