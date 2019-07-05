struct DiscretePath
    ticks::Vector{Float64}
    ticklabels::Vector{String}
    positions::Vector{Float64}
    points::Matrix{Float64}
end

# Initialization
function DiscretePath(kDict0::PointDict, named_path::Vector{String}; num_points::Int=60)

    ticklabels = [kDict0.label[key] for key in named_path]

    ticks, positions, points = names_to_path(named_path, kDict0.coord, num_points)

    DiscretePath(ticks, ticklabels, positions, points)
end

function DiscretePath(kDict0::PointDict; num_points::Int=60)

    named_path = kDict0.default_path
    DiscretePath(kDict0, named_path; num_points=num_points)
end

# Point iterator
function eachpoint(kPoints::DiscretePath)
    eachcol(kPoints.points)
end

# Manipulators
function scaled_ticks(kPoints::DiscretePath; start=0.0, length=1.0)
    k_ticks = kPoints.ticks
    start .+ k_ticks/maximum(k_ticks)  * (length-start)
end

##################################################################################
##################################################################################
##################################################################################

function exportdata(filename, ks::DiscretePath)
    h5open(filename, "w") do file
        write(file, "ticks", ks.ticks)
        write(file, "ticklabels", ks.ticklabels)
        write(file, "positions", ks.positions)
        write(file, "points", ks.points)
    end
end

##################################################################################
##################################################################################
##################################################################################

function names_to_coords(names, point_dictionary)
    [point_dictionary[name] for name in names]
end


function coords_to_points(B, coords)
    [Vector(B * c) for c in coords]
end


function names_to_path(names, point_dictionary, N; B=[1.0 0.0; 0.0 1.0])
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

    return tick_coordinates, parametric_path, hcat(path_array...)
end
