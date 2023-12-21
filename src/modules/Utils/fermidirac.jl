heaviside(x) = ifelse(real.(x) < 0, zero(x), ifelse(real.(x) > 0, one(x), oftype(x, 0.5)))

fermidirac_T(ϵ; T::Number=0.01) = 1.0 ./ (exp.(ϵ ./ T) .+ 1)
fermidirac_0(ϵ) = heaviside.(-ϵ)
fermidirac(ϵ; T = 0, μ = 0) = (T > 0) ? fermidirac_T(ϵ .- μ; T = T) : fermidirac_0(ϵ .- μ)

function groundstate_sumk(ϵs_k::AbstractVector, μ::Number=0; T::Number=0)
    sum(fermidirac(ϵs_k .- μ; T=T) .* ϵs_k)
end