dim(A::AbstractMatrix, x) = size(A,1)
dim(f::Function, x::Number) = size(f(x), 1)
dim(f::Function, x::AbstractVector) = size(f(first(x)), 1)
dim(f::Function, x::AbstractMatrix) = size(f(first(eachcol(x))), 1)
dim(h::Hops, x) = size(first(values(h)),1)


function expvalf(𝑶::AbstractMatrix)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶 * ψ)
    f
end

function expvalf(𝑶::Function)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)
    f
end

function expvalf(𝑶::Hops)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)
    f
end