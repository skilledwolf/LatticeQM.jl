mutable struct LabeledPoints
    coord::Dict{String, AbstractVector{Float64}}
    label::Dict{String, String}
    defaultpath::Vector{String}
end

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
