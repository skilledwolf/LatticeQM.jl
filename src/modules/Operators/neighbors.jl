

import LinearAlgebra: norm
import ..Utils: @scalar2vector
import ..TightBinding: Hops, addspin

function getneighborhops(lat::Lattice; cellrange::Int=3)

    δt(r) = (0-DEFAULT_PRECISION<r<float(cellrange)/sqrt(3)) ? r : -1
    @scalar2vector δt

    diffhops = Dict(δR=>round.(real(M); digits=9) for (δR,M) in Hops(lat, δt; cellrange=cellrange))

    diff = filter(x->x>-0.5, sort(unique(vcat((unique(M) for (δR,M) in diffhops)...))))
    
    diff, [Hops(δR=>complex(float(M.==D)) for (δR,M) in diffhops if norm(M.==D)>0) for D in diff]

end


function getneighborhops(lat, t0, t1=0, t2=0; cellrange=4, spin=true)

    dists, neighborhops = getneighborhops(lat; cellrange=cellrange)
    hops1, hops2, hops3 = neighborhops[1:3]

    hops = t0 * hops1 

    if t1 != 0
        hops += t1 * hops2
    end
    if t2 != 0
        hops += t2 * hops3
    end

    if spin
        hops = addspin(hops, :spinhalf)
    end

    Hops(R=>complex(M) for (R,M) in hops)
end