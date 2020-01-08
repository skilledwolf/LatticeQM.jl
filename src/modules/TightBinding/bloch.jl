
BlochPhase(k,δL) = exp(1.0im * 2 * π * complex(dot(k,δL)))

fouriersum(hoppings, k::AbstractVector) = sum(complex(t) .* BlochPhase(k, δL) for (δL,t) in hoppings)

get_bloch(hoppings::AnyHops) = k -> fouriersum(hoppings, k)
