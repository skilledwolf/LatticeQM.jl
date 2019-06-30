using Base

####################################################################################################
## K-points, paths and discrete paths
####################################################################################################

struct KPoints{T}
    positions::ElasticArray{T,2}
    names::Vector{String}
    labels::Vector{String}
end
KPoints{T}(dim::Int64=1) where {T<:Real} = KPoints{T}(Array{T}(undef,dim,0), [], [])
KPoints(args...) = KPoints{Float64}(args...)

struct KPath
    names::Vector{String}
    defnum::Int64
    default::Bool
end
KPath(names; defnum=100, default=true) = KPath(names, defnum, default)

struct DiscreteKPath{T<: Real}
    kticks::Vector{String}
    kparametric::Vector{T}
    kpoints::Matrix{T}
end

####################################################################################################
####################################################################################################
####################################################################################################

function Base.:+(kpoints1::KPoints, kpoints2::KPoints)
    append!(kpoints1.positions, kpoints2.positions)
    append!(kpoints1.names, kpoints2.names)
    append!(kpoints1.labels, kpoints2.labels)

    kpoints1
end

####################################################################################################
####################################################################################################
####################################################################################################


function names_to_coords(names, point_dictionary)
    [point_dictionary[name] for name in names]
end


function coords_to_points(B, coords)
    [Vector(B * c) for c in coords]
end


function names_to_path(names, point_dictionary, B, N)
    construct_path(coords_to_points(B, names_to_coords(names, point_dictionary)), N)
end


function construct_path(points, N)

    num_points = length(points)
    differences = [points[i]-points[i-1] for i=2:num_points]
    differences_len = [norm(diff) for diff in differences]
    differences = differences ./ differences_len
    cumtotal = cumsum(differences_len)
    pushfirst!(cumtotal, 0.0)

    path_array = []

#     collect(range(0; length=N, stop=cumtotal[end]))
    parametric_path = collect(range(0; length=N, stop=cumtotal[end]))

    for s in range(0; length=N, stop=cumtotal[end])
        del_s = s .- cumtotal

        for i=2:num_points
            if del_s[i-1] >= 0 && del_s[i] <= 0
                push!(path_array, points[i-1] + (del_s[i-1] .* differences[i-1]))
                continue
            end
        end
    end

    tick_coordinates = cumtotal


    return tick_coordinates, parametric_path, path_array
end