const DEFAULT_PRECISION = sqrt(eps())::Float64

import ..Structure
import ..Structure.Lattices
import ..Structure.Lattices: Lattice

import ..TightBinding

import ..TightBinding: autoconversion

TightBinding.Hops(lat::Lattice, args...; kwargs...) = gethops(lat, args...; kwargs...)
TightBinding.addhops!(hops::Hops, lat::Lattice, t::Function, args...; kwargs...) = TightBinding.addhops!(hops, gethops(lat, t, args...; kwargs...))

"""
    gethops(lat::Lattice, t::Function; cellrange=1, format=:auto, vectorized=false)

Iterates over pairs of orbitals/atom positions (r1,r2) in lattice `lat` and evaluates
the hopping elements t(r1+R,r2) for each lattice vector R.

By default, `vectorized=false`. For huge systems use `vectorized=true` and make 
sure the hopping function t accepts matrices as arguments.
The keyword argument `format` can be `:dense` or `:sparse`. For `:auto`, small systems 
will be dense and huge problems are assumed to be sparse.

Returns the hopping elements in the format
`Dict(R => t_R)`

"""
function gethops(lat::Lattice, args...; format=:auto, kwargs...)
    hops = autoconversion(Hops(), format)
    hops!(hops, lat, args...; kwargs...)
end

import ..Structure.Lattices: getneighbordict

function hops!(hops::Hops, lat::Lattice, t::Function; cellrange=2, vectorized=false, kwargs...)
    R = Lattices.allpositions(lat)
    neighbor_dict = getneighbordict(lat, cellrange)
    if vectorized
        vectorizedhops!(hops, R, neighbor_dict, t; kwargs...)
    else
        hops!(hops, R, neighbor_dict, t; kwargs...)
    end
    TightBinding.trim!(hops)
    hops
end

###############################################################################
# Main routines for gethops(...)
###############################################################################

asserthopdim(t0::Number) = 1
asserthopdim(t0::AbstractMatrix) = size(t0,1)

function vectorizedhops!(hops::Hops, R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function) #, format=:auto    
    R2 = similar(R)
    for (δL,δa) in neighbors
        R2 .= R .+ δa
        hops[δL] = t(R2, R)
        hops[-δL] = deepcopy(hops[δL]') # create the Hermitian conjugates
    end
    hops
end

import ..TightBinding: zero_matrix

function hops!(hops::Hops{K,T}, R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; kwargs...) where {K,T}
    N = size(R,2)
    d = asserthopdim(t(R[:,1]))::Int
    V  = Matrix{ComplexF64}(undef, (d, d)) # preallocate memory for the hopping matrix

    for (δL,δa) in neighbors
        # Initialize a zero matrix of type T
        hops[δL] = zero_matrix(typeof(hops), d*N, d*N)
        hoppingmatrix!(hops[δL], V, R .+ δa, R, t; kwargs...) # heavy lifting
        hops[-δL] = deepcopy(hops[δL]') # create the Hermitian conjugates
    end
    hops
end

function hoppingmatrix!(M::AbstractMatrix{ComplexF64}, 
                             V::Array{ComplexF64}, 
                             Ri::Matrix{Float64}, 
                             Rj::Matrix{Float64}, 
                             t::Function; 
                             precision::Float64=DEFAULT_PRECISION,
                             maxsize::Int=0)
    
    d = size(V, 2) # bond dimension
    N = size(Ri, 2) # number of atoms
    element_count = 0

    for i = 1:N, j = 1:N
        @views V[1:d, 1:d] .= t(Ri[:,i], Rj[:,j])  # Update V directly
        for i0 = 1:d, j0 = 1:d
            val = V[i0, j0]
            if abs(val) >= precision
                row = (i - 1) * d + i0
                col = (j - 1) * d + j0
                M[row, col] = val  # Update the sparse matrix directly
                element_count += 1
            end

            @assert maxsize < 1  || element_count <= maxsize "The number of elements in the sparse matrix exceeded the specified maxsize."
        end
    end
    S
end