CappedYukawa(r::AbstractVector{Float64}; kwargs...) = CappedYukawa(norm(r); kwargs...)
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
Hubbard(r::AbstractVector{Float64}; kwargs...) = Hubbard(norm(r); kwargs...)
Hubbard(r::Float64; a=0.9, U=1.0) = U * heaviside(a-r)

build_CappedYukawa(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->CappedYukawa(r; kwargs...); mode=mode, format=format)
build_Hubbard(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->Hubbard(r; kwargs...); mode=mode, format=format)
