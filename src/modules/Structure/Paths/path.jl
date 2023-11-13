using LaTeXStrings: latexstring

using OrderedCollections: OrderedDict

struct DiscretePath
    points::Matrix{Float64}
    ticks::OrderedDict{Float64,String}
    positions::Vector{Float64} # linear positions along path
end

# tickindices(ks::DiscretePath) = keys(ks.ticks)
ticks(ks::DiscretePath) = collect(keys(ks.ticks))
ticklabels(ks::DiscretePath) = collect(values(ks.ticks))

function Base.show(io::IO, ks::DiscretePath)
    println(io, "Discrete Path: ", join(latexstring.(ticklabels(ks)), "→"), "  (",size(ks.points,2)," points)")
end

# Initialization
DiscretePath(ks::Matrix, ticks::Vector, ticklabels::Vector{String}, pos::Vector) = DiscretePath(ks, OrderedDict(f=>s for (f,s)=zip(ticks,ticklabels)),pos)

function DiscretePath(kdict::LabeledPoints, named_path::Vector{String}; B::AbstractMatrix=1.0I, num_points::Int=60)
    
    ticklabels, points = kdict(named_path)
    points, ticks, positions = getpath(B * points; num=num_points, B=B)

    DiscretePath(points, ticks, ticklabels, positions)
end

function DiscretePath(kdict::LabeledPoints; kwargs...)

    named_path = kdict.defaultpath
    DiscretePath(kdict, named_path; kwargs...)
end

# Interface for iteration and item access
Base.eachcol(ks::DiscretePath) = Base.eachcol(ks.points)
Base.size(ks::DiscretePath, args...) = Base.size(ks.points, args...)

Base.iterate(ks::DiscretePath) = length(ks)>0 ? (first(eachcol(ks.points)),Base.firstindex(ks)) : nothing
Base.iterate(ks::DiscretePath, state) = (state+=1; state>lastindex(ks) ? nothing : (ks[state],state))
Base.length(ks::DiscretePath) = Base.size(ks,2)
# Base.eltype(ks::DiscretePath) = typeof(ks[firstindex(ks)])

Base.values(ks::DiscretePath) = Base.eachcol(ks.points)
Base.firstindex(ks::DiscretePath) = first(CartesianIndices(ks.points))[2]
Base.lastindex(ks::DiscretePath) = last(CartesianIndices(ks.points))[2]
Base.get(ks::DiscretePath, args...) = Base.get(ks.points, args...)
Base.getindex(ks::DiscretePath, args...) = Base.getindex(ks.points, args...)
Base.setindex!(ks::DiscretePath,v,args...) = Base.setindex!(ks.points,v,args...)

################################################################################
################################################################################


# Point iterator
# Define proper iterators for each input type
# points(ks::DiscretePath) = ks.points
# points(ks::T) where {T<:AbstractMatrix} = ks
# points(ks::T2) where {T2<:AbstractVector{<:AbstractVector}} = hcat(ks...)


# function sumk(f_k::Function, ks)
#     sum(f_k(k) for k=eachcol(ks))/size(ks,2)
# end

# import ..regulargrid
# sumk(f_k::Function; klin::Int) = sumk(f_k, regulargrid(;nk=klin^2))


################################################################################
################################################################################

# Manipulators
function scaleticks(kPoints::DiscretePath; start=0.0, length=1.0)
    k_ticks = ticks(kPoints)
    start .+ k_ticks/maximum(k_ticks)  * (length-start)
end


################################################################################
################################################################################

import DelimitedFiles

DelimitedFiles.writedlm(ks::DiscretePath) = DelimitedFiles.writedlm(".", ks)

function DelimitedFiles.writedlm(path::String, ks::DiscretePath)
    @assert ispath(path) "Error: '$path' is not a valid path."
    writedlm(dirname(path)*"/ticks.out", ticks(ks))
    writedlm(dirname(path)*"/ticklabels.out", ticklabels(ks))
    writedlm(dirname(path)*"/positions.out", ks.positions)
    writedlm(dirname(path)*"/points.out", ks.points)
end

