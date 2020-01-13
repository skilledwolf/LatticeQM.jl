heaviside(x::Real) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

fermidirac_T(ϵ::Real; T::Real=0.01) = 1.0/(exp(ϵ/T)+1)
fermidirac_0(ϵ::Real) = heaviside(-ϵ)

function fermidirac(ϵ::Real; T::Real=0.01)
    if T==0.0
        x = fermidirac_0(ϵ)
    else
        x = fermidirac_T(ϵ; T=T)
    end

    x
end
