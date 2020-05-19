"""
    blockmatrix(mat, I, J, V)
    blockmatrix!(mat, I, J, V)

Writes a sparse matrix from block matrices.

N (Int): Size of the resulting N-times-N square matrix
I,J (Vector{Int}): Row/Column coordinates
vec_of_mats (Vector{Matrix}): Vector of block matrices (must be equal sized in the current implementation)
"""
function blockmatrix!(mat::AbstractMatrix, I::Vector{Int}, J::Vector{Int}, vec_of_mats::AbstractVector{<:AbstractMatrix})
    @assert length(I) == length(J) == length(vec_of_mats)
    @assert size(unique([I J]; dims=1),1)== length(I) "Block coordinates must be unique."
    
    for (i,j,m) in zip(I,J,vec_of_mats)
        blockmatrix!(mat, i, j, m)
    end

    mat
end

function blockmatrix!(mat::AbstractMatrix, i::Int, j::Int, m::AbstractMatrix)

    s = size(m)
    X0 = CartesianIndex(s[1]*(i-1),s[2]*(j-1))

    for X in CartesianIndices(m)
        mat[X0+X] = m[X]
    end

    mat
end

function blockmatrix(N::Int, I::Vector{Int}, J::Vector{Int}, vec_of_mats::AbstractVector{<:AbstractMatrix})

    mat = spzeros(N,N)
    blockmatrix!(mat, I, J, vec_of_mats)

    mat
end

####################################################################################################
####################################################################################################
####################################################################################################

using ..Structure.Lattices: supercellpoints

"""
    superlattice(hops, periods)

Turn a given hopping model into a superlattice model by copying cells and hoppings.
This method can be useful as preparation before adding modulations that change the
periodicity of the model.
"""
superlattice(hops::AnyHops, v::Vector{Int}, args...; kwargs...) = superlattice(hops, Matrix(Diagonal(v)), args...; kwargs...)
superlattice(hops::AnyHops, M::Matrix{Int}) = superlattice(hops, M, (r,R)->1)
function superlattice(hops::AnyHops, M::Matrix{Int}, phasefunc::Function) where {T<:Number}

    coordinates = supercellpoints(M)
    
    count = size(coordinates, 2)
    D = count  * hopdim(hops) # size of superlattice hopping matrices

    sneighbors = hcat([[i;j] for i=-2:2 for j=-2:2]...)

    shops = Hops(L => spzeros(ComplexF64, D,D) for L=eachcol(sneighbors))

    basislookup = Dict(L => i for (L, i) in zip(eachcol(coordinates), 1:count))

    for (j,a) in enumerate(eachcol(coordinates))
        for (δa, t) in hops
            v = round.(a + δa)
            n = findfirst(i->haskey(basislookup, v + M' * sneighbors[:,i]), 1:size(sneighbors,2))

            L = (n==nothing) ? (print(eltype(coordinates)); error("Error: target atom at coordinate $v not found.")) : sneighbors[:,n] # throw error if the required position does not exist
            i = basislookup[v+ M'*L]

            z = phasefunc([i,j], δa) # phase

            blockmatrix!(shops[L], i, j, t .* z) # write to the corresponding block of the corresponding superlattice hopping matrix
        end
    end

    shops
end