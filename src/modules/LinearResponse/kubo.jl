
import ..Utils: fermidirac


# function kuboterm(ω::Number, e1::Float64, e2::Float64, M::Number; T::Float64=0.01, μ::Float64=0) #todo: test
#     (fermidirac(e1; T=T, μ=μ)-fermidirac(e2, T=T, μ=μ)) * M / (ω + (e1-e2)) # + 1im*Γ
# end

import LinearAlgebra

function kubo!(output::Vector{ComplexF64}, Ω::AbstractVector, ϵs::AbstractVector, 
    states::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix;
    μ::Float64=0.0, Γ::Float64=0.025, T::Float64=0.1)

    M = (states' * A * states) .* transpose(states' * B' * states)
    fd = fermidirac(ϵs; T=T, μ=μ)

    N = length(ϵs)
    it = ((i,j) for i=1:N, j=1:N if i!=j)

    for (i_,j_)=it
        output .+= -1im .* (fd[i_]-fd[j_]) .* M[i_,j_] ./ (Ω .+ 1im*Γ .+ (ϵs[i_]-ϵs[j_]))
    end

    output
end

kubo(Ω::Number, args...; kwargs...) = kubo([Ω], args...; kwargs...)
function kubo(Ω::AbstractVector, args...; kwargs...)
    output = complex(zero(Vector(Ω)))
    kubo!(output, Ω, args...; kwargs... )
end