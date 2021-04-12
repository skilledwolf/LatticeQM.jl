using LinearAlgebra: norm, I, transpose

getpath(points::AbstractMatrix; kwargs...) = getpath(collect(eachcol(points)); kwargs...)
function getpath(points::AbstractVector; num::Int=100, B=I)

    num_points = length(points)

    if num_points==1 # if there's only one point there is no use in discretizing anything ;-)
        return transpose(B * inv(transpose(B) * B)) * hcat(points...), [0.0], [0.0]
    end
    
    differences = [points[i]-points[i-1] for i=2:num_points]
    differences_len = [norm(diff) for diff in differences]
    differences = differences ./ differences_len
    cumtotal = cumsum(differences_len)
    pushfirst!(cumtotal, 0.0)

    path_array = []
    parametric_path = collect(range(0; length=num, stop=cumtotal[end]))
    for s in range(0; length=num, stop=cumtotal[end])
        del_s = s .- cumtotal

        for i=2:num_points
            if del_s[i-1] >= 0 && del_s[i] <= 0
                push!(path_array, points[i-1] + (del_s[i-1] .* differences[i-1]))
                continue
            end
        end
    end

    tick_coordinates = cumtotal
    return transpose(B * inv(transpose(B) * B)) * hcat(path_array...), tick_coordinates, parametric_path
end