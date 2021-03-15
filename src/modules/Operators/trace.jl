
import ..TightBinding: Hops, AnyHops, zerokey

import LinearAlgebra: tr, Diagonal

trace(A::AbstractMatrix) = tr(A)

# Structure.Paths.sumk(getbloch(ρ_sol), kgrid)
trace(hops::Hops) = tr(hops[zerokey(hops)]) #sum(tr(m) for m in values(hops))
trace(h1::Hops,h2::Hops) = sum(tr(h1[L]*h2[L]) for L=intersect(keys(h1),keys(h2)))

trace(hops::Hops, A::AbstractMatrix) = tr(hops[zerokey(hops)]*A)

trace(hops::Hops, lat, name::String) = trace(hops, getoperator(lat, name))
trace(hops::Hops, lat, names::AbstractVector{String}) = [trace(hops, o) for o in getoperator(lat, names)]

trace(hops::Hops, operators::AbstractVector) = [trace(hops,o) for o=operators]

function density(ρ::Hops)
    d = trace(ρ)
    @assert isapprox(imag(d),0; atol=sqrt(eps())) "Complex-valued density? Seems unphysical, please check."
    real(d)
end


# This is the correct one to use for contraction with "density matrix"
expval(h1::Hops) = sum(sum(h1[L]) for L=keys(h1))
expval(h1::Hops, h2::Hops) = sum(sum(h1[L].*h2[L]) for  L=intersect(keys(h1),keys(h2)))
expval(h1::Hops, h2::AbstractMatrix) = expval(h1, Hops(zerokey(h1)=>h2))
expval(h1::AbstractMatrix, h2::Hops) = expval(Hops(h1), h2)
expval(hops::Hops, name::String, lat::Lattice) = expval(hops, getoperator(lat, name)')
expval(hops::Hops, operators::AbstractVector, args...) = [expval(hops,o, args...) for o=operators]

magnetization(ρ, lat::Lattice) = expval(ρ, ["sx", "sy", "sz"], lat)
magnetization(ρ, A::AbstractMatrix, lat::Lattice) = expval(ρ, [A*s for s=getoperator(lat, ["sx", "sy", "sz"])])
magnetization(ρ, A::AbstractVector{<:AbstractMatrix}, lat::Lattice) = [magnetization(ρ, a, lat) for a=A]
magnetization(ρ, A::AbstractVector{String}, lat::Lattice) = [magnetization(ρ, getoperator(lat,a), lat) for a=A]


import ..Structure
import ..Structure.Lattices: Lattice

function localdensity(ρ, lat::Lattice)

    N = Structure.countorbitals(lat)
    D = zeros(ComplexF64, N)

    P(i) = kron(Diagonal([float(i==j) for j=1:N]), Diagonal(ones(2)))

    for i=1:N
        D[i] = expval(ρ, P(i))
    end

    D
end

function localmagnetization(ρ, lat::Lattice)

    N = Structure.countorbitals(lat)
    M = zeros(ComplexF64, 3,N)

    P(i) = kron(Diagonal([float(i==j) for j=1:N]), Diagonal(ones(2)))

    for i=1:N
        M[:,i] = real.(magnetization(ρ, P(i), lat))
    end

    M
end