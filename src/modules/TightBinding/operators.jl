

import ..Spectrum: expvalf, dim 

Spectrum.dim(h::Hops, x=0) = size(first(values(h)),1)

function Spectrum.expvalf(𝑶::Hops)
    f(k, ψ, ϵ) = real.(ψ' * 𝑶(k) * ψ)
    f
end