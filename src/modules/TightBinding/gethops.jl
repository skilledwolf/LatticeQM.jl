###############################################################################
# Wrapper for custom types to gethops(...)
###############################################################################

import ..Structure
import ..Structure.Lattices: Lattice

addhops!(hops::Hops, lat::Lattice, t::Function; kwargs...) = addhops!(hops, gethops(lat, t; kwargs...))

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
function gethops(lat::Lattice, t::Function; cellrange=1, format=:auto, precision::Float64=sqrt(eps()), kwargs...)# where {T<:AbstractMatrix{Float64}}
    # Get neighbor cells
    neighbors = Structure.getneighborcells(lat, cellrange; halfspace=true, innerpoints=true, excludeorigin=false)
    # Iterate the hopping function over orbital pairs and neighbors
    gethops(lat, neighbors, t; precision=precision, format=format, kwargs...)
end

import ..Utils: padvec

function gethops(lat::Lattice, neighbors::Vector{Vector{Int}}, t::Function; kwargs...)
    R = Structure.allpositions(lat)
    d = size(R,1)

    A = Structure.getA(lat)
    neighbor_dict = Dict(δL => padvec(A*δL,d) for δL in neighbors)

    gethops(R, neighbor_dict, t; kwargs...)
end


###############################################################################
# Main routines for gethops(...)
###############################################################################

asserthopdim(t0::Number) = 1
asserthopdim(t0::AbstractMatrix) = size(t0,1)

function gethops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; vectorized=false, format=:auto, kwargs...)

    if vectorized # indicates that the function t accepcts the call signature t(R1::Matrix,R2::Matrix)
        return getvectorizedhops(R, neighbors, t; format=format, kwargs...)
    end

    if format==:auto
        format = decidetype(size(R,1))
    end

    if format==:dense
        getdensehops(R, neighbors, t; kwargs...)
    elseif format==:sparse
        getsparsehops(R, neighbors, t; kwargs...)
    else
        error("Format `$format` does not exist. Choose `:auto`, `:dense` or `sparse`.")
    end
end

using Distributed

function getvectorizedhops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; precision::Float64=1e-6, format=:auto)
    N = size(R,2)
    hops = Hops()

    for (δL,δa) in neighbors
        Ri = R.+δa
        hops[δL] = t(Ri, R)
        hops[-δL] = hops[δL]' # create the Hermitian conjugates
    end

    ensuretype(hops, format)
end

function getdensehops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; precision::Float64=1e-6)
    N = size(R,2)
    d = asserthopdim(t(R[:,1]))::Int

    # Preallocate memory: avoid unnecessary allocations
    M  = Matrix{ComplexF64}(undef, (d*N, d*N))

    hops = Hops()
    for (δL,δa) in neighbors
        densehoppingmatrix!(M, R.+δa, R, t)
        hops[δL] = copy(M)
        hops[-δL] = hops[δL]' # create the Hermitian conjugates
    end

    hops
end

function getsparsehops(R::Matrix{Float64}, neighbors::Dict{Vector{Int},Vector{Float64}}, t::Function; precision::Float64=1e-6)
    N = size(R,2)
    d = asserthopdim(t(R[:,1]))::Int

    # Preallocate memory: important for huge sparse matrices
    maxind = (N>MAX_DENSE) ? round(Int, MAX_DIAGS * N) : N^2 # semi-arbitrary limit for dense allocation
    IS = Vector{Int}(undef, maxind*d^2)
    JS = similar(IS)
    VS = similar(IS, ComplexF64)
    V  = Matrix{ComplexF64}(undef, (d, d))

    hops = Hops()
    for (δL,δa) in neighbors
        hops[δL] = sparsehoppingmatrix!(IS, JS, VS, V, R.+δa, R, t; precision=precision) # heavy lifting
        hops[-δL] = hops[δL]' # create the Hermitian conjugates
    end

    hops
end

function densehoppingmatrix!(M::Array{ComplexF64}, Ri::Matrix{Float64}, Rj::Matrix{Float64}, t::Function)
    N = size(Ri,2) # number of atoms
    d = div(size(M, 2), N)
    @assert mod(size(M, 2), d)==0 "Incompatible dimensions."

    # Iterate over atom pairs and calculate the (possibly matrix-valued) hopping amplitudes
    for i=1:N, j=1:N
        I = (i-1)*d; J = (j-1)*d
        @views M[I+1:I+d, J+1:J+d] .= t(Ri[:,i], Rj[:,j])
    end

    M
end


import SparseArrays: sparse

function sparsehoppingmatrix!(IS::Vector{Int}, JS::Vector{Int}, VS::Array{ComplexF64}, V::Array{ComplexF64}, Ri::Matrix{Float64}, Rj::Matrix{Float64}, t::Function; precision::Float64)

    d = size(V, 2) # bond dimension
    N = size(Ri,2) # number of atoms
    maxind = div(length(IS),d^2) # preallocated memory

    count = 0 # counter for added matrix elements

    for i=1:N,j=1:N
        @views V[1:d, 1:d] .= t(Ri[:,i], Rj[:,j])

        for i0=1:d, j0=1:d
            if abs(V[i0, j0]) < precision
                continue
            end

            count = count+1

            if count > maxind 
                error("The tight-binding matrix is not sparse enough for this constructor.")
            end
            
            IS[count] = (i-1)*d + i0
            JS[count] = (j-1)*d + j0
            VS[count] = V[i0, j0]
        end
    end

    sparse(IS[1:count], JS[1:count], VS[1:count], N*d, N*d)
end
