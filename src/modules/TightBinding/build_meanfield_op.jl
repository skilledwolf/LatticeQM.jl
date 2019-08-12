# Scalar generators
CappedYukawa(r::AbstractVector{Float64}; kwargs...) = CappedYukawa(norm(r); kwargs...)
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
Hubbard(r::AbstractVector{Float64}; kwargs...) = Hubbard(norm(r); kwargs...)
Hubbard(r::Float64; a=0.5, U=1.0) = U * heaviside(a-r)


# Lattice operators
function get_Hubbard(lat, neighbors=[[0;0]]; mode=:nospin, format=:auto, kwargs...)
    """
    returns Dict(δL => Matrix(V(r_i-r_j+δL))_ij) where δL are vectors that
    connect unit cells. The set of δL's (in units of lattice vectors) is specified by 'neighbors'.
    """
    ee_exchange = get_hops(lat, neighbors, r->Hubbard(r; kwargs...))
    extend_space!(ee_exchange, mode)

    ee_exchange
end

function get_CappedYukawa(lat, neighbors=[[i;j] for i=-1:1 for j=-1:1]; mode=:nospin, kwargs...)
    ee_exchange = get_hops(lat, neighbors, r->CappedYukawa(r; kwargs...))
    extend_space!(ee_exchange, mode)

    ee_exchange
end

build_CappedYukawa(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->CappedYukawa(r; kwargs...); mode=mode, format=format)
build_Hubbard(lat; mode=:nospin, format=:auto, kwargs...) = build_H(lat, r->Hubbard(r; kwargs...); mode=mode, format=format)
