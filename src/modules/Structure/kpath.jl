function names_to_coords(names, point_dictionary)
    [point_dictionary[name] for name in names]
end


function coords_to_points(B, coords)
    [Vector(B * c) for c in coords]
end


function names_to_path(names, point_dictionary, N; B=I)
    construct_path(coords_to_points(B, names_to_coords(names, point_dictionary)), N; B=B)
end


function construct_path(points, N; B=I)

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

    return tick_coordinates, parametric_path, inv(B) * hcat(path_array...)
end
