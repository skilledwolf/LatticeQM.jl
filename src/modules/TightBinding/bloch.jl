fourierphase(k,δL) = exp(1.0im * 2 * π * complex(dot(k,δL)))

fouriersum(hoppings, k::AbstractVector) = sum(complex(t) .* fourierphase(k, δL) for (δL,t) in hoppings)

getbloch(hoppings) = k -> fouriersum(hoppings, k)

function fouriersum(hoppings, k::Real, d::Int)
    N=length(zerokey(hoppings))

    K = unique(collect(key[(1:N).!=d] for key=keys(hoppings)))
    newhops = Dict(key => zero(getzero(hoppings)) for key in K)

    for (key,h)=hoppings
        newhops[key[(1:N).!=d]] .+= complex(h) .* fourierphase([k],[key[d]])
    end

    newhops
end

getbloch(hoppings,d::Int) = k -> fouriersum(hoppings,k,d)


###################################################################################################
# Legacy maps
###################################################################################################
@legacyalias getbloch get_bloch
@legacyalias fourierphase BlochPhase
