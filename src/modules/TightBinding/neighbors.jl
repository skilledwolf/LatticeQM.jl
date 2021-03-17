

import ..Utils: @scalar2vector

function getneighborhops(lat::Lattice; cellrange::Int=3)

    δt(r) = (0-1e-10<r<float(cellrange)/sqrt(3)) ? r : -1
    @scalar2vector δt

    diffhops = Dict(δR=>round.(real(M); digits=9) for (δR,M) in gethops(lat, δt; cellrange=cellrange))

    diff = filter(x->x>-0.5, sort(unique(vcat([unique(M) for (δR,M) in diffhops]...))))
    
    diff, [Hops(δR=>complex(float(M.==D)) for (δR,M) in diffhops if norm(M.==D)>0) for D in diff]

end