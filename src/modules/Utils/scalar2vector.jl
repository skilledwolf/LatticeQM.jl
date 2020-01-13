macro scalar2vector(f0, N=3)
"""
This is macro is a wrapper that takes as input a function f0(x::Float64) and adds a new dispatch
f0(r1::Vector, r2::Vector) = f0(norm(r1-r2)) while making sure that r1 and r2 do not exceed length N.
"""
    return quote
        function $(esc(f0))(r1::AbstractVector{Float64}, r2::Float64=0.0; kwargs...)
            n = min(length(r1), $N)
            $(esc(f0))(norm(r1[1:n].-r2); kwargs...)
        end
        function $(esc(f0))(r1::AbstractVector{Float64}, r2::AbstractVector{Float64}; kwargs...)
            n = min(length(r1), $N)
            $(esc(f0))(norm(r1[1:n].-r2[1:n]); kwargs...)
        end

        $(esc(f0))
    end
end