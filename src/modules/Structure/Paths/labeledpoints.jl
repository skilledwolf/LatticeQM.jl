mutable struct LabeledPoints
    points::Dict{String,Vector{Float64}}
    labels::Dict{String, String}
    defaultpath::Vector{String}
end

labels(lpoints::LabeledPoints, S::AbstractVector{String}) = collect(lpoints.labels[s] for s=S)
points(lpoints::LabeledPoints, S::AbstractVector{String}) = hcat((lpoints.points[s] for s=S)...)

(lpoints::LabeledPoints)(s::String) = lpoints[s]
(lpoints::LabeledPoints)(S::AbstractVector{String}) = (labels(lpoints, S), points(lpoints, S))#?

# Interface for iteration and item access
import Base
# Base.get(lpoints::LabeledPoints, args...) = (Base.get(lpoints.labels, args...), Base.get(lpoints.points, args...))
Base.haskey(lpoints::LabeledPoints, args...) = Base.haskey(lpoints.points, args...)
Base.length(lpoints::LabeledPoints) = Base.length(lpoints.points)
# Base.eltype(lpoints::LabeledPoints) = Base.eltype(lpoints.points)

Base.values(lpoints::LabeledPoints) = Base.values(lpoints.points)
Base.keys(lpoints::LabeledPoints) = Base.keys(lpoints.points)
Base.getindex(lpoints::LabeledPoints, args...) = Base.getindex(lpoints.points, args...)
Base.setindex!(lpoints::LabeledPoints, v, args...) = Base.setindex!(lpoints.points, v, args...)


# Initialization
function LabeledPoints(names::Vector{String}, coord::Vector{T}, label::Vector{String}, default::Vector{String}) where {T<:AbstractVector{Float64}}
    coord_dict = Dict(a=>b for (a,b) in zip(names, coord))
    labels_dict = Dict(a=>b for (a,b) in zip(names, label))

    LabeledPoints(coord_dict, labels_dict, default)
end

function LabeledPoints(names::Vector{String}, coord::AbstractMatrix{Float64}, label::Vector{String}, default::Vector{String})
    LabeledPoints(names, Vector(collect(eachcol(coord))), label, default)
end

LabeledPoints(d::Int=0) = LabeledPoints(Vector{String}(), zeros(d,d), Vector{String}(), Vector{String}())


# Modification
function Base.:append!(kDict0::LabeledPoints, name::String, coord::AbstractVector{Float64}, label::String)
    kDict0.coord[name] = coord
    kDict0.label[name] = label

    return nothing
end

# Display
function Base.:display(kDict0::LabeledPoints)
    display(kDict0.coord)
    display(kDict0.label)
    display(kDict0.defaultpath)
end
