
using ..Utils: fermidirac

function kubo(ω, e1::Float64,e2::Float64, M::Float64; T::Float64=0.01) #todo: test
    (fermidirac(e1; T=T)-fermidirac(e2, T=T))/(e1-e2) * M ./ (ω - (e1-e2)) # + 1im*Γ
end

function kubo(Ω, energies::AbstractVector, states::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix;
    μ::Float64=0.0, Γ::Float64=0.025, T::Float64=0.1) # todo: test

    result = complex(zero(Ω))

    M = (states' * A * states) .* (states * B * states') # todo: test!

    Ωshift = Ω + 1im*Γ # finite lifetime shift / pole shift

    for (i_,e1) = enumerate(energies.-μ), (j_,e2) = enumerate(energies.-μ)
        if i_ != j_
                result += kubo(Ωshift, e1, e2, M[i_,j_]; kwargs...)
        end
    end

    result
end