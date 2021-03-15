struct DiscretePath
    ticks::Vector{Float64}
    ticklabels::Vector{String}
    positions::Vector{Float64} # linear positions along path
    points::Matrix{Float64}
end

using LaTeXStrings: latexstring

function Base.show(io::IO, ks::DiscretePath)
    println(io, "Discrete Path: ", join(latexstring.(ks.ticklabels), "→"), "  (",size(ks.points,2)," points)")
end

using HDF5
using ...DummySave
function DummySave.save!(file, ks::DiscretePath) # file should be hdf5 output stream
    g = g_create(file, "path")

    g["ticks"] = ks.ticks
    g["ticklabels"] = ks.ticklabels
    g["positions"] = ks.positions
    g["points"] = ks.points
end

DummySave.save(data::DiscretePath, filename::String="kpath.h5") = DummySave.save_wrapper(data, filename)

################################################################################
################################################################################

# Initialization
function DiscretePath(kDict0::LabeledPoints, named_path::Vector{String}; B::AbstractMatrix=1.0I, num_points::Int=60)

    ticklabels = [kDict0.label[key] for key in named_path]

    ticks, positions, points = names2path(named_path, kDict0.coord, num_points; B=B)

    DiscretePath(ticks, ticklabels, positions, points)
end

function DiscretePath(kDict0::LabeledPoints; kwargs...)

    named_path = kDict0.defaultpath
    DiscretePath(kDict0, named_path; kwargs...)
end

################################################################################
################################################################################

import ...Utils

# Point iterator
# Define proper iterators for each input type
const kIterable = Union{DiscretePath, <:AbstractMatrix{Float64}, <:AbstractVector{T1}} where {T1<:AbstractVector{Float64}}
eachpoint(kPoints::DiscretePath) = eachcol(kPoints.points)
eachpoint(ks::T) where {T<:AbstractMatrix{Float64}} = eachcol(ks)
eachpoint(ks::T2) where {T1<:AbstractVector{Float64},T2<:AbstractVector{T1}} = ks
points(kPoints::DiscretePath) = kPoints.points
points(ks::T) where {T<:AbstractMatrix{Float64}} = ks
points(ks::T2) where {T2<:AbstractVector{<:AbstractVector{Float64}}} = hcat(ks...)

function sumk(f_k::Function, ks::kIterable)
    ks = points(ks)
    sum(f_k(k) for k=eachcol(ks))/size(ks,2)
end

function sumk(f_k::Function; klin::Int)
    ks = Utils.regulargrid(;nk=klin^2)
    sum(f_k(k) for k=eachcol(ks))/size(ks,2)
end

################################################################################
################################################################################

# Manipulators
function scaleticks(kPoints::DiscretePath; start=0.0, length=1.0)
    k_ticks = kPoints.ticks
    start .+ k_ticks/maximum(k_ticks)  * (length-start)
end
