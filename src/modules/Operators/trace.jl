
import ..TightBinding: Hops, AbstractHops, zerokey

import LinearAlgebra: tr, Diagonal

trace(A::AbstractMatrix) = tr(A)

# Structure.Paths.sumk(getbloch(ρ_sol), kgrid)
trace(hops::Hops) = tr(hops[zerokey(hops)]) #sum(tr(m) for m in values(hops))

# (1/Nk) Σ_k tr[A(k)B(k)] = Σ_L tr[A_L B_{−L}] in the Bloch convention
# M(k) = Σ_L M_L e^{2πik·L}. The L block of A must contract against the −L
# block of B; pairing L with L instead computes Σ_L tr[A_L B_L], which is
# wrong for every inter-cell block and flips the sign of transpose-odd
# onsite operators (e.g. sy).
trace(h1::Hops,h2::Hops) = sum(_blocktrace(h1[L], h2[-L]) for L=keys(h1) if haskey(h2,-L))

# tr(A*B) without forming the product.
_blocktrace(A, B) = sum(transpose(B) .* A)

trace(hops::Hops, A::AbstractMatrix) = tr(hops[zerokey(hops)]*A)

trace(hops::Hops, lat, name::String) = trace(hops, getoperator(lat, name))
trace(hops::Hops, lat, names::AbstractVector{String}) = [trace(hops, o) for o in getoperator(lat, names)]

trace(hops::Hops, operators::AbstractVector) = [trace(hops,o) for o=operators]

function density(ρ::Hops)
    d = trace(ρ)
    @assert isapprox(imag(d),0; atol=sqrt(eps())) "Complex-valued density? Seems unphysical, please check."
    real(d)/hopdim(ρ)
end


# This is the correct one to use for contraction with "density matrix":
# ⟨O⟩ = (1/Nk) Σ_k tr[ρ(k)O(k)] = Σ_L tr[ρ_L O_{−L}]  (see `trace` above).
expval(h1::Hops) = sum(sum(h1[L]) for L=keys(h1))
expval(h1::Hops, h2::Hops) = trace(h1, h2)
expval(h1::Hops, h2::AbstractMatrix) = expval(h1, Hops(zerokey(h1)=>complex(h2)))
expval(h1::AbstractMatrix, h2::Hops) = expval(Hops(complex(h1)), h2)
expval(hops::Hops, name::String, lat::Lattice) = expval(hops, getoperator(lat, name))
expval(hops::Hops, operators::AbstractVector, args...) = [expval(hops,o, args...) for o=operators]

magnetization(ρ, lat::Lattice) = expval(ρ, ["sx", "sy", "sz"], lat)
magnetization(ρ, A::AbstractMatrix, lat::Lattice) = expval(ρ, [A*s for s=getoperator(lat, ["sx", "sy", "sz"])])
magnetization(ρ, A::AbstractVector{<:AbstractMatrix}, lat::Lattice) = [magnetization(ρ, a, lat) for a=A]
magnetization(ρ, A::AbstractVector{String}, lat::Lattice) = [magnetization(ρ, getoperator(lat,a), lat) for a=A]


import ..Structure
import ..Structure.Lattices
import ..Structure.Lattices: Lattice

"""Site densities, summed over the internal states belonging to each orbital."""
function localdensity(ρ::Hops, lat::Lattice)
    N = Lattices.countorbitals(lat)
    block = ρ[zerokey(ρ)]
    size(block, 1) == size(block, 2) && size(block, 1) % N == 0 ||
        throw(DimensionMismatch("Density matrix must contain an equal number of internal states per orbital"))
    n = size(block, 1) ÷ N
    result = zeros(ComplexF64, N)
    for i in 1:N, a in 1:n
        j = n * (i - 1) + a
        result[i] += block[j, j]
    end
    result
end

"""Site magnetizations in the spin-half basis (spin is the fast index)."""
localmagnetization(ρ::Hops, lat::Lattice) = localobservables(ρ, lat)[2:4, :]

"""Return site density and Pauli x/y/z expectations as rows, in O(N) work."""
function localobservables(ρ::Hops, lat::Lattice)
    N = Lattices.countorbitals(lat)
    block = ρ[zerokey(ρ)]
    size(block) == (2N, 2N) || throw(DimensionMismatch("Spin observables require two states per orbital"))
    result = Matrix{ComplexF64}(undef, 4, N)
    for i in 1:N
        up, down = 2i - 1, 2i
        a, b = block[up, up], block[up, down]
        c, d = block[down, up], block[down, down]
        result[1, i] = a + d
        result[2, i] = b + c
        result[3, i] = im * (b - c)
        result[4, i] = a - d
    end
    result
end
