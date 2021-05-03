import ..TightBinding
import ..TightBinding: Hops
import ..Structure.Lattices

"""
    Returns current operators [J_α,...] for each spatial coordinate α=1,...,D
    given a generating function t for the hopping elements of the respective Hamiltonian.
"""
function getcurrentoperators(lat::Lattice, t::Function; mode=:nospin, kwargs...) # todo: test
    D = Lattices.latticedim(lat)
    ts = [(r1::AbstractVector,r2::AbstractVector -> (-1im)*(r1-r2)[i]* t(r1,r2; kwargs...)) for i=1:D]

    [Hops(lat, ts[i]; kwargs...) for i=1:D]
end

"""
    Returns current operators [J_α,...] for each spatial coordinate α=1,...,D
    given a Hamiltoinan of hopping elements `hops` and lattice structure `la`.

    This should be the preferred way of generating a current operator for a given hopping Hamiltonian.
"""
function getcurrentoperators(lat::Lattice, hops::AbstractMatrix) # todo: test
    N = Lattices.countorbitals(lat)
    d = 2
    D = size(hops,1)

    @assert D>=N && mod(D,N) == 0
    n = div(D,N) # for spin-1/2 we would have n=2

    r = Lattices.allpositions(lat)[1:2,:]

    currentoperators = [deepcopy(hops) for l_=1:d]
    for l_=1:d
        for i_=1:N, j_=1:N # todo: test!
            currentoperators[l_][1+(i_-1)*n:n+(i_-1)*n,1+(j_-1)*n:n+(j_-1)*n] .*= -1im * (r[:,i_]-r[:,j_])[l_] * ones(n,n)
        end
    end

    currentoperators
end




"""
    Returns current operators [J_α,...] for each spatial coordinate α=1,...,D
    given a Hamiltoinan of hopping elements `hops` and lattice structure `la`.

    This should be the preferred way of generating a current operator for a given hopping Hamiltonian.
"""
function getcurrentoperators(lat::Lattice, hops::Hops) # todo: test
    N = Lattices.countorbitals(lat)
    # d = latticedim(lat)
    d = Lattices.spacedim(lat)
    D = TightBinding.hopdim(hops)

    @assert D>=N && mod(D,N) == 0
    n = div(D,N) # for spin-1/2 we would have n=2

    r = Lattices.positions(lat)

    Ls = Dict(L => Lattices.getA(lat) * L for L in keys(hops))

    currentoperators = [deepcopy(hops) for l_=1:d]

    for l_=1:d
        for (L,A) in Ls
            r0 = r.+A
            currentoperators!(currentoperators[l_][L], r0[l_,:], r[l_,:])
        end
    end

    currentoperators
end

function currentoperators!(M::AbstractMatrix, r1, r2)

    if issparse(M)
        # M[:] = sparse(M)[:]
        return currentoperators_sparse!(M, r1, r2)
    else 
        return currentoperators_dense!(M, r1, r2)
    end
end

function currentoperators_dense!(M::AbstractMatrix, r1, r2)
    N = length(r1)
    D = size(M,2)
    n = div(D,N)

    for j_=1:N 
        δR = r1.-r2[j_]
        for i_=1:N
            M[1+(i_-1)*n:n+(i_-1)*n,1+(j_-1)*n:n+(j_-1)*n] .*= -1im * δR[i_]
        end
    end
    M
end

function currentoperators_sparse!(M::SparseMatrixCSC, r1, r2)
    N = length(r1)
    D = size(M,2)
    n = div(D,N)

    rows = rowvals(M)
    # vals = nonzeros(M)

    myind(i) = round(Int,1+(i-1-rem(i-1,n))/n)

    for j_ in 1:size(M, 2)
        δR = r1.-r2[myind(j_)]
        for r in nzrange(M, j_)
            i_= rows[r] 
            M[i_,j_] *= -1im * δR[myind(i_)]
            # v_ = vals[r]
        end
    end
end
