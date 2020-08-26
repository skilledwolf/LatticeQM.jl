
using ..TightBinding: zerokey

trace(A::AbstractMatrix) = tr(A)

# Structure.Paths.sumk(getbloch(ρ_sol), kgrid)
trace(hops::AnyHops) = tr(hops[zerokey(hops)]) #sum(tr(m) for m in values(hops))
trace(h1::AnyHops,h2::AnyHops) = sum(tr(h1[L]*h2[L]) for L=intersect(keys(h1),keys(h2)))

trace(hops::AnyHops, A::AbstractMatrix) = tr(hops[zerokey(hops)]*A)

trace(hops::AnyHops, lat, name::String) = trace(hops, getoperator(lat, name))
trace(hops::AnyHops, lat, names::AbstractVector{String}) = [trace(hops, o) for o in getoperator(lat, names)]

trace(hops::AnyHops, operators::AbstractVector) = [trace(hops,o) for o=operators]

function density(ρ::AnyHops)
    d = trace(ρ)
    @assert isapprox(imag(d),0; atol=sqrt(eps())) "Complex-valued density? Seems unphysical, please check."
    real(d)
end

# This is the correct one to use for contraction with "density matrix"
expval(h1::AnyHops) = sum(sum(h1[L]) for L=keys(h1))
expval(h1::AnyHops, h2::AnyHops) = sum(sum(h1[L].*h2[L]) for  L=intersect(keys(h1),keys(h2)))
expval(h1::AnyHops, h2::AbstractMatrix) = expval(h1, Hops(h2))
expval(h1::AbstractMatrix, h2::AnyHops) = expval(Hops(h1), h2)
expval(hops::AnyHops, name::String, lat::Lattice) = expval(hops, getoperator(lat, name)')
expval(hops::AnyHops, operators::AbstractVector, args...) = [expval(hops,o, args...) for o=operators]

magnetization(ρ, lat::Lattice) = expval(ρ, ["sx", "sy", "sz"], lat)
magnetization(ρ, A::AbstractMatrix, lat::Lattice) = expval(ρ, [A*s for s=getoperator(lat, ["sx", "sy", "sz"])])
magnetization(ρ, A::AbstractVector{<:AbstractMatrix}, lat::Lattice) = [magnetization(ρ, a, lat) for a=A]
magnetization(ρ, A::AbstractVector{String}, lat::Lattice) = [magnetization(ρ, getoperator(lat,a), lat) for a=A]

using ..Structure: countorbitals

function localmagnetization(ρ, lat::Lattice)

    N = countorbitals(lat)
    M = zeros(3,N)

    P(i) = kron(Diagonal([float(i==j) for j=1:N]), Diagonal(ones(2)))

    for i=1:N
        M[:,i] = real.(magnetization(ρ, P(i), lat))
    end

    M
end