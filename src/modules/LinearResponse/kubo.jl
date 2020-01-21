
using ..Utils: fermidirac

function kubo(ω, e1::Float64, e2::Float64, M::Number; T::Float64=0.01) #todo: test
    (fermidirac(e1; T=T)-fermidirac(e2, T=T))/(e1-e2) * M ./ (ω .- (e1 - e2)) # + 1im*Γ
end

function kubo(Ω, ϵs::AbstractVector, states::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix;
    μ::Float64=0.0, Γ::Float64=0.025, T::Float64=0.1) # todo: test

    @assert maximum(abs.(imag(ϵs)))<1e-15   "Complex energies might be unphysical!"
    ϵs = real(ϵs)

    result = complex(zero(Ω))

    M = (states' * A * states) .* transpose(states' * B' * states) # todo: test!
#     M = [dot(conj(ψ1),A*ψ2) * dot(conj(ψ2),B'*ψ1) for ψ1=eachcol(states), ψ2=eachcol(states)]

    Ωshift = Ω .+ 1im*Γ # finite lifetime shift / pole shift

    for (i_,e1) = enumerate(ϵs.-μ), (j_,e2) = enumerate(ϵs.-μ)
        if i_ != j_
                result .+= kubo(Ωshift, e1, e2, M[i_,j_]; T=T)
        end
    end

    -1im * result
end