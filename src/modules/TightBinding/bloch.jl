fourierphase(k,δL) = exp(1.0im * 2 * π * complex(dot(k,δL)))

fouriersum(hoppings, k::AbstractVector) = sum(complex(t) .* fourierphase(k, δL) for (δL,t) in hoppings)

getbloch(hoppings) = k -> fouriersum(hoppings, k)


###################################################################################################
# Legacy maps
###################################################################################################
@legacyalias getbloch get_bloch
@legacyalias fourierphase BlochPhase
