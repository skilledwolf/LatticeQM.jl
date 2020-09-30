using ..TightBinding: gethops, addspin
using ..Utils: heaviside, @scalar2vector

# Functions with scalar arguments
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))
Hubbard(r::Float64; a=0.5, U=1.0) = U * heaviside(a-r)

# Functions with vector arguments
@scalar2vector CappedYukawa
@scalar2vector Hubbard

# Lattice operators
function gethubbard(lat, neighbors=[[0;0]]; mode=:nospin, format=:auto, kwargs...)
    """
    returns Dict(δL => Matrix(V(r_i-r_j+δL))_ij) where δL are vectors that
    connect unit cells. The set of δL's (in units of lattice vectors) is specified by 'neighbors'.
    """
    t(args...) = Hubbard(args...; kwargs...)
    ee_exchange = gethops(lat, neighbors, t; format=format)

    addspin(ee_exchange, mode)
end

# function getcappedyukawa(lat; cellrange=1, mode=:nospin, format=:auto, kwargs...)
#     t(args...) = CappedYukawa(args...; kwargs...)
#     ee_exchange = gethops(lat, neighbors, t; cellrange=cellrange, format=format)

#     addspin(ee_exchange, mode)
# end

function getcappedyukawa(lat, args...; mode=:nospin, k0=1.0, U=1.0, kwargs...)
    t(args0...) = CappedYukawa(args0...; k0=k0, U=U)
    ee_exchange = gethops(lat, args..., t; kwargs...)

    addspin(ee_exchange, mode)
end

# build_CappedYukawa(lat; mode=:nospin, format=:auto, kwargs...) = gethamiltonian(lat, r->CappedYukawa(r; kwargs...); mode=mode, format=format)
# gethamiltonianubbard(lat; mode=:nospin, format=:auto, kwargs...) = gethamiltonian(lat, r->Hubbard(r; kwargs...); mode=mode, format=format)



###################################################################################################
# Backwards compatibility
###################################################################################################
export get_Hubbard
@legacyalias gethubbard get_Hubbard

export get_Hubbard
@legacyalias getcappedyukawa get_CappedYukawa
