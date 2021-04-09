
function expvalf(𝑶::AbstractMatrix)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶 * ψ)
    f
end

function expvalf(𝑶::Function)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)
    f
end