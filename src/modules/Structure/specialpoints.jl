mutable struct PointDict
    coord::Dict{String, AbstractVector{Float64}}
    label::Dict{String, String}
    default_path::Vector{String}
end

# Initialization
# PointDict() = PointDict(Dict{String,AbstractVector{Float64}}(), Dict{String,String}(), Vector{String}([]))

function PointDict(names::Vector{String}, coord::Vector{T}, label::Vector{String}, default::Vector{String}) where {T<:AbstractVector{Float64}}
    coord_dict = Dict(a=>b for (a,b) in zip(names, coord))
    labels_dict = Dict(a=>b for (a,b) in zip(names, label))

    PointDict(coord_dict, labels_dict, default)
end

function PointDict(names::Vector{String}, coord::AbstractMatrix{Float64}, label::Vector{String}, default::Vector{String})
    PointDict(names, Vector(collect(eachcol(coord))), label, default)
end

# Modification
function Base.:append!(kDict0::PointDict, name::String, coord::AbstractVector{Float64}, label::String)
    kDict0.coord[name] = coord
    kDict0.label[name] = label

    return nothing
end

# Display
function Base.:display(kDict0::PointDict)
    display(kDict0.coord)
    display(kDict0.label)
    display(kDict0.default_path)
end
