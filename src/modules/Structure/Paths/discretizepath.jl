function names2coordinates(names, dict)
    [dict[name] for name in names]
end

function coordinates2points(B, coords)
    [Vector(B * c) for c in coords]
end

function names2path(names, dict, N; B=I)
    path(coordinates2points(B, names2coordinates(names, dict)), N; B=B)
end

function path(points::AbstractVector, N; B=I)

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

    return tick_coordinates, parametric_path,  transpose(B * inv(transpose(B) * B)) * hcat(path_array...)
end
