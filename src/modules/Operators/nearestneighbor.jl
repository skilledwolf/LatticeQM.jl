import ..TightBinding: Hops, addhops!
import LinearAlgebra: norm

"""
    DistanceWindowHopping(a, halfwidth, t0)

Callable hopping rule: returns `t0` when the bond length is in
`(a - halfwidth, a + halfwidth)`, otherwise `zero(t0)`.

The struct is callable with three signatures (matching `gethops`'s
expectations):
    h(r::Real)                                   # scalar distance
    h(r1::AbstractVector, r2::Real=0.0)          # uses norm(r1) (truncated to 3D)
    h(r1::AbstractVector, r2::AbstractVector)    # uses norm(r1 - r2) (truncated to 3D)

`t0` may be any number type; the return type is `T = typeof(t0)` regardless of
the branch, so downstream `gethops`/`hops!` see a type-stable function. Every
call to `nearestneighbor!` produces a value of the *same* concrete struct type
(parameterised by `T`), so the `gethops` → `hops!` → `hoppingmatrix!` pipeline
compiles once per `T` instead of once per closure.
"""
struct DistanceWindowHopping{T<:Number}
    a::Float64
    halfwidth::Float64
    t0::T
end

# Truncation to 3D matches the previous @scalar2vector(t, 3) behaviour:
# lattices may carry extra (non-spatial) coordinates (e.g. "layer") in
# components 4+, which must not contribute to bond length.
@inline function _spatial_norm(r1::AbstractVector{Float64}, r2::Float64)
    n = min(length(r1), 3)
    norm(@view(r1[1:n]) .- r2)        # broadcast scalar offset, matches old macro
end
@inline function _spatial_norm(r1::AbstractVector{Float64}, r2::AbstractVector{Float64})
    n = min(length(r1), length(r2), 3)
    norm(@view(r1[1:n]) .- @view(r2[1:n]))
end

(h::DistanceWindowHopping{T})(r::Real) where {T} =
    (h.a + h.halfwidth > r > h.a - h.halfwidth) ? h.t0 : zero(T)

(h::DistanceWindowHopping{T})(r1::AbstractVector{Float64}, r2::Float64=0.0) where {T} =
    h(_spatial_norm(r1, r2))

(h::DistanceWindowHopping{T})(r1::AbstractVector{Float64}, r2::AbstractVector{Float64}) where {T} =
    h(_spatial_norm(r1, r2))

getnearestneighborhops(args...; kwargs...) = nearestneighbor!(Hops(), args...; kwargs...)

"""
    nearestneighbor!(hops, lat, t0=-1.0; a=1.0, halfwidth=0.01, kwargs...)

Add nearest-neighbour hopping `t0` (between sites at distance ≈ `a`) to `hops`.

Backed by `DistanceWindowHopping(a, halfwidth, t0)`; pass that struct directly
to `addhops!` if you need the same hopping rule across multiple lattices
without redundant compilation.
"""
function nearestneighbor!(hops, lat, t0=-1.0; a=1.0, halfwidth=0.01, kwargs...)
    addhops!(hops, lat, DistanceWindowHopping(float(a), float(halfwidth), t0); kwargs...)
end
