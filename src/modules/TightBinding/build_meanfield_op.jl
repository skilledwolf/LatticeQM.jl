CappedYukawa(r::AbstractVector{Float64}; kwargs...) = CappedYukawa(norm(r); kwargs...)
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))

build_CappedYukawa(lat; mode=:nospin, kwargs...) = build_H(lat, r->CappedYukawa(r; kwargs...); mode=mode)
