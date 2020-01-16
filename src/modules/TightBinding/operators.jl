dim(A::AbstractMatrix, x) = size(A,1)
dim(f::Function, x::Number) = size(f(x), 1)
dim(f::Function, x::AbstractVector) = size(f(first(x)), 1)
dim(f::Function, x::AbstractMatrix) = size(f(first(eachcol(x))), 1)
dim(h::AnyHops, x) = size(first(values(h)),1)


function expvalf(𝑶::AbstractMatrix)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶 * ψ)
    f
end

expvalf(𝑶::AnyHops) = expvalf(getbloch(𝑶))

function expvalf(𝑶::Function)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)
    f
end