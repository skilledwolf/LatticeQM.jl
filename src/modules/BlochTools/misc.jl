heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

fermidirac_T(ϵ::AbstractFloat; T::AbstractFloat=0.01) = 1.0/(exp(ϵ/T)+1)
fermidirac_0(ϵ::AbstractFloat) = heaviside(-ϵ)

function fermidirac(ϵ::AbstractFloat; T::AbstractFloat=0.01)
    if T==0.0
        x = fermidirac_0(ϵ)
    else
        x = fermidirac_T(ϵ; T=T)
    end

    x
end
