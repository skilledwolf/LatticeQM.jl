import ..TightBinding
import ..TightBinding: Hops
import ..Structure.Lattices

# Current operators are J_α(k) = ∂H(k)/∂k_α, with α running over the lattice
# directions only. We use `latticedim`, not `spacedim` — for a 2D lattice
# embedded in 3D space `spacedim` would produce a third operator that is
# identically zero (the Bloch phase has no k_z component to differentiate).

"""
    getcurrentoperators(lat, hops::Hops)

Return the current operators `[J_1, ..., J_d]` (one per lattice direction)
for the tight-binding Hamiltonian represented by `hops` on lattice `lat`.

Each `J_α` is an `Hops` whose `(i,j)` element at offset `R` is
`-i (r_i + A·R - r_j)_α · h_{ij}(R)`, i.e. ∂/∂k_α of the Bloch Hamiltonian
with the standard position-included gauge.
"""
function getcurrentoperators(lat::Lattice, hops::Hops)
    N = Lattices.countorbitals(lat)
    d = Lattices.latticedim(lat)
    D = TightBinding.hopdim(hops)

    @assert D >= N && mod(D, N) == 0
    n = div(D, N)   # spin or other internal multiplicity

    r = Lattices.positions(lat)
    Ls = Dict(L => Lattices.getA(lat) * L for L in keys(hops))

    Js = [deepcopy(hops) for _ in 1:d]
    for α in 1:d, (L, A) in Ls
        r0 = r .+ A
        _apply_position_factor!(Js[α][L], r0[α, :], r[α, :])
    end
    Js
end

"""
    getcurrentoperators(lat, t::Function; kwargs...)

Build current operators directly from a hopping function `t(r1, r2; ...)`.
Returns one `Hops` per lattice direction. `kwargs` are forwarded both to `t`
and to the `Hops(lat, ...)` constructor (e.g. `cellrange`, `format`).
"""
function getcurrentoperators(lat::Lattice, t::Function; kwargs...)
    d = Lattices.latticedim(lat)
    [Hops(lat, (r1, r2) -> -1im * (r1 - r2)[α] * t(r1, r2; kwargs...);
          kwargs...) for α in 1:d]
end


# ---------------------------------------------------------------------------
# Internal: in-place position-factor application
# ---------------------------------------------------------------------------
# Multiplies each (i,j) block of `M` by `-i · (r1[i] - r2[j])`, the
# k-derivative of the Bloch phase factor in the position-included gauge.
# Block size is `n × n` per orbital pair, where `n = size(M, 1) ÷ length(r1)`
# (1 for spinless, 2 for spin-1/2, etc.).

_apply_position_factor!(M::AbstractMatrix, r1, r2) =
    SparseArrays.issparse(M) ? _apply_position_factor_sparse!(M, r1, r2) :
                               _apply_position_factor_dense!(M, r1, r2)

function _apply_position_factor_dense!(M::AbstractMatrix, r1, r2)
    N = length(r1)
    n = div(size(M, 2), N)
    @inbounds for j in 1:N
        δR = r1 .- r2[j]
        for i in 1:N
            @views M[1+(i-1)*n : i*n, 1+(j-1)*n : j*n] .*= -1im * δR[i]
        end
    end
    M
end

function _apply_position_factor_sparse!(M::SparseMatrixCSC, r1, r2)
    N = length(r1)
    n = div(size(M, 2), N)
    blockindex(i) = div(i - 1, n) + 1   # 1-based orbital index from row/col

    rows = rowvals(M)
    @inbounds for col in 1:size(M, 2)
        δR = r1 .- r2[blockindex(col)]
        for r in nzrange(M, col)
            i_block = blockindex(rows[r])
            M[rows[r], col] *= -1im * δR[i_block]
        end
    end
    M
end
