heaviside(x) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

fermidirac_T(ϵ; T=0.01) = 1.0/(exp(ϵ/T)+1)
fermidirac_0(ϵ) = heaviside(-ϵ)

fermidirac(ϵ::Number; T=0) = (T>0) ? fermidirac_T(ϵ; T=T) : fermidirac_0(ϵ)

fermidirac(ϵs::AbstractVector; kwargs...) = [fermidirac(ϵ) for ϵ=ϵs]

