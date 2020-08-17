@legacyalias getnearestneighborhops get_NN_hops
getnearestneighborhops(args...; kwargs...) = nearestneighbor!(Hops(), args...; kwargs...)

function nearestneighbor!(hops, lat, t0=1.0; a=1.0, selector=(:))
    t(r) = (a+0.01>r[:]>a-0.01) ? t0 : 0.0
    @scalar2vector t

    addhops!(hops, lat, t)
end