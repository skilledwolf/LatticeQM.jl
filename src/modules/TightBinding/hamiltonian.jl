gethamiltonian(args...; mode=:nospin, kwargs...) = getbloch(gethops(args...; kwargs...); mode=mode)

###############################################################################
# Wrapper for custom types to gethops(...)
###############################################################################

using ..Structure: getneighborcells

addhops!(hops::AnyHops, lat::Lattice, t::Function; kwargs...) = addhops!(hops, gethops(lat, t; kwargs...))

function gethops(lat::Lattice, t::Function; cellrange=1, format=:auto, precision::Float64=1e-6, kwargs...)# where {T<:AbstractMatrix{Float64}}

    # Get lattice neighbors
#     neighbors = [[i;j] for i=-1:1 for j=-1:1 if i+j>=0 && i>=0]
    neighbors = getneighborcells(lat, cellrange; halfspace=true, innerpoints=true, excludeorigin=false)

    gethops(lat, neighbors, t; precision=precision, format=format, kwargs...)
end

using ..Utils: padvec

function gethops(lat::Lattice, neighbors::Vector{Vector{Int}}, t::Function; kwargs...)
    R = allpositions(lat)
    d = size(R,1)

    A = getA(lat)
    neighbor_dict = Dict(δL => padvec(A*δL,d) for δL in neighbors)

    gethops(R, neighbor_dict, t; kwargs...)
end


###############################################################################
# Main routines for gethops(...)
###############################################################################

const MAX_DENSE = 500
const MAX_DIAGS = 100

function gethops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; vectorized=false, kwargs...)

    ####
    # The keyword specialized indicates that the function t accepcts the signature t(R1,R2) where
    # R1 and R2 are matrices. This can be crucial for performance when dealing with huge unit cells
    # (such as for twister bilayer graphene).
    # I might make this the default, but then it should be well-documented.
    ####

    if vectorized
        return gethops_vectorized(R, neighbors, t; kwargs...)
    end

    ####
    # Switch between matrix-valued t and scalar t functions.
    # For large unit cells matrix-valued t is a bad idea!
    ###
    v0 = t(R[:,1])
    if isa(v0, Number)
        d = 1
    elseif isa(v0, AbstractMatrix)
        d = size(v0,1)
    else
        error("Hopping function t must return number or matrix of number")
    end

    return gethops(R, neighbors, t, d; kwargs...)
end

using Distributed

function gethops_vectorized(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; precision::Float64=1e-6, format=:auto)

    N = size(R,2)
    hops = Hops()

    for (δL,δa) in neighbors
        Ri = R.+δa
        hops[δL] = t(Ri, R)
        hops[-δL] = hops[δL]' # create the Hermitian conjugates
    end

    decidetype(hops, format)
end

function gethops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function, d::Int; precision::Float64=1e-6, format=:auto)
    N = size(R,2)

    # Preallocate memory: important for huge sparse matrices
    maxind = (N>MAX_DENSE) ? round(Int, MAX_DIAGS * N) : N^2 # semi-arbitrary limit for dense allocation
    IS = Vector{Int}(undef, maxind*d^2)
    JS = similar(IS)
    VS = similar(IS, ComplexF64)
    V  = Matrix{ComplexF64}(undef, (d, d))

    hops = Hops()
    for (δL,δa) in neighbors
        hops[δL] = hoppingmatrix!(IS, JS, VS, V, R.+δa, R, t; precision=precision) # heavy lifting
        hops[-δL] = hops[δL]' # create the Hermitian conjugates
    end

    decidetype(hops, format)
end

function hoppingmatrix!(IS::Vector{Int}, JS::Vector{Int}, VS::Array{ComplexF64}, V::Array{ComplexF64}, Ri::Matrix{Float64}, Rj::Matrix{Float64}, t::Function; precision::Float64)

    d = size(V, 2) # bond dimension
    N = size(Ri,2) # number of atoms
    maxind = div(length(IS),d^2) # preallocated memory

    count = 0 # counter for added matrix elements

    @fastmath for i=1:N,j=1:N
        @views V[1:d, 1:d] .= t(Ri[:,i], Rj[:,j])

        for i0=1:d, j0=1:d
            if abs(V[i0, j0]) < precision
                continue
            end
            count = count+1
            IS[count] = (i-1)*d + i0
            JS[count] = (j-1)*d + j0
            VS[count] = V[i0, j0]
        end
    end

    sparse(IS[1:count], JS[1:count], VS[1:count], N*d, N*d)
end


###################################################################################################
# Backwards compatibility
###################################################################################################
@legacyalias gethops get_hops