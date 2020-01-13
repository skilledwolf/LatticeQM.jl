mutable struct MeanFieldData{T<:AbstractArray{<:Number,2}}
    ρ::Dict{Vector{Int},T}
    ϵGS::Union{Nothing, Array{Float64,3}}
    converged::Bool
    error::Float64
end
