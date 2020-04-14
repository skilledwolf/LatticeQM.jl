
"""
    Returns current operators [J_α,...] for each spatial coordinate α=1,...,D
    given a generating function t for the hopping elements of the respective Hamiltonian.
"""
function getcurrentoperators(lat::Lattice, t::Function; mode=:nospin, kwargs...) # todo: test
    D = latticedim(lat)
    ts = [(r1::AbstractVector,r2::AbstractVector -> (-1im)*(r1-r2)[i]* t(r1,r2; kwargs...)) for i=1:D]

    [gethops(lat, ts[i]; kwargs...) for i=1:D]
end

"""
    Returns current operators [J_α,...] for each spatial coordinate α=1,...,D
    given a Hamiltoinan of hopping elements `hops` and lattice structure `la`.

    This should be the preferred way of generating a current operator for a given hopping Hamiltonian.
"""
function getcurrentoperators(lat::Lattice, hops::AnyHops) # todo: test
    N = countorbitals(lat)
    d = latticedim(lat)
    D = hopdim(hops)

    @assert D>=N && mod(D,N) == 0
    n = div(D,N) # for spin-1/2 we would have n=2

    r = positions(lat)

    Ls = Dict(L => getA(lat) * L for L in keys(hops))

    currentoperators = [deepcopy(hops) for l_=1:d]
    for (L,A) in Ls, l_=1:d
        for i_=1:N, j_=1:N # todo: test!
            currentoperators[l_][L][1+(i_-1)*n:n+(i_-1)*n,1+(j_-1)*n:n+(j_-1)*n] .*= -1im * (r[:,i_].+A-r[:,j_])[l_] * ones(n,n)
        end
    end

    currentoperators
end

"""
    Returns current operators [J_α,...] for each spatial coordinate α=1,...,D
    given a Hamiltoinan of hopping elements `hops` and lattice structure `la`.

    This should be the preferred way of generating a current operator for a given hopping Hamiltonian.
"""
function getcurrentoperators(lat::Lattice, hops::AbstractMatrix) # todo: test
    N = countorbitals(lat)
    d = 2
    D = size(hops,1)

    @assert D>=N && mod(D,N) == 0
    n = div(D,N) # for spin-1/2 we would have n=2

    r = allpositions(lat)[1:2,:]

    currentoperators = [deepcopy(hops) for l_=1:d]
    for l_=1:d
        for i_=1:N, j_=1:N # todo: test!
            currentoperators[l_][1+(i_-1)*n:n+(i_-1)*n,1+(j_-1)*n:n+(j_-1)*n] .*= -1im * (r[:,i_]-r[:,j_])[l_] * ones(n,n)
        end
    end

    currentoperators
end