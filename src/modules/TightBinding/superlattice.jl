function block_sparse_quad(N::Int, I::Vector{Int}, J::Vector{Int}, vec_of_mats::AbstractVector{T}) where {T<:AbstractMatrix}
"""
    Builds a sparse matrix by from block matrices.

    N (Int): Size of the resulting N-times-N square matrix
    I,J (Vector{Int}): Row/Column coordinates
    vec_of_mats (Vector{Matrix}): Vector block matrices.

    Note that all vectors I,J,vec_of_mats must have the same length.
"""
    mat = spzeros(N,N)

    for (i,j,m) in zip(I,J,vec_of_mats)
        s = size(m)
        I0 = CartesianIndex(s[1]*(i-1),s[2]*(j-1))

        for I in CartesianIndices(m)
            mat[I0+I] = m[I]
        end
    end

    mat
end

####################################################################################################
####################################################################################################
####################################################################################################

function superlattice_hops(hops::Dict{Vector{Int},AbstractMatrix{T}},lattice_vecs::Matrix{Int}) where {T<:Number}
"""
    Turn a given hopping model into a superlattice model by trivially copying cells.
    This method can be useful as preparation before adding modulations that change the
    periodicity of the model.

"""

    intlattice = points_within_supercell(lattice_vecs')
    basis_lookup = Dict{Vector{Int}, Int}(vec => ind for (vec, ind) in zip(eachcol(intlattice), 1:(size(intlattice)[2])))

    orbdim = size(collect(values(hops))[1])[2]
    uc_size = size(intlattice)[2]
    N = uc_size * orbdim

    cellbasis = 1:uc_size

    supercellneighbors = lattice_vecs * hcat([[i;j] for i=-1:1 for j=-1:1]...)

    suplat_jump_IJ   = Dict{Vector{Int}, ElasticArray{Int,2}}()
    suplat_jump_VALS = Dict{Vector{Int}, ElasticArray{Matrix{Complex{Float64}}}}()

    for suplatvec in eachcol(supercellneighbors)
        suplat_jump_IJ[suplatvec] = ElasticArray{Int}(undef, 2, 0)
        suplat_jump_VALS[suplatvec] = ElasticArray{Matrix{Complex{Float64}}}(undef, 0)
    end

    for j=cellbasis
        latpos = intlattice[:,j]
        for (lathop, t) in hops
            found = false
            for suplathop in eachcol(supercellneighbors)
                try
                    i = basis_lookup[latpos + lathop + suplathop]
                    append!(suplat_jump_IJ[suplathop], [i; j])
                    append!(suplat_jump_VALS[suplathop], [t])

                    found = true
                    break

                    catch error
                    if isa(error, KeyError)
                        continue
                    else
                        throw(error)
                    end
                end
            end
            if !found
                error("target atom not found: "*string(latpos+lathop))
            end
        end
    end

    suplat_hoppings = Dict(
        key=>block_sparse_quad(N,suplat_jump_IJ[key][1,:],suplat_jump_IJ[key][2,:],suplat_jump_VALS[key])
        for key in keys(suplat_jump_IJ)
    )
end