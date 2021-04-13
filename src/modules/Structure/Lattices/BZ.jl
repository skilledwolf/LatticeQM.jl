import ..Paths: getpath, DiscretePath

mutable struct BrillouinZone
    B::AbstractMatrix{Float64}
    points::Dict{String,Vector{Float64}}
end

(BZ::BrillouinZone)(s::String) = BZ[s]
(BZ::BrillouinZone)(S::AbstractVector{String}) = hcat(collect(BZ[s] for s=S)...)#?
function (BZ::BrillouinZone)(S::AbstractVector{String}, num::Int)
    ks, ticks, pos = getpath(BZ.B*BZ(S); num=num, B=BZ.B)

    DiscretePath(ks,pos,ticks,S)
end

# Interface for iteration and item access
import Base
Base.get(BZ::BrillouinZone, args...) = Base.get(BZ.points, args...)
Base.haskey(BZ::BrillouinZone, args...) = Base.haskey(BZ.points, args...)
Base.length(BZ::BrillouinZone) = Base.length(BZ.points)
Base.eltype(BZ::BrillouinZone) = Base.eltype(BZ.points)

Base.values(BZ::BrillouinZone) = Base.values(BZ.points)
Base.keys(BZ::BrillouinZone) = Base.keys(BZ.points)
Base.setindex!(BZ::BrillouinZone, v, args...) = Base.setindex!(BZ.points, v, args...)

