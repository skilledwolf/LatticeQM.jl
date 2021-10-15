using ..TightBinding: Hops, addspin
using ..Utils: heaviside, @scalar2vector

using ..TightBinding: zerokey

# Functions with scalar arguments
CappedYukawa(r::Float64; k0=1.0, U=1.0) = U/(k0*r*exp(k0*r)+exp(k0*r))
Hubbard(r::Float64; a=0.5, U=1.0) = U * heaviside(a-r)
Tanh(r::Float64; a=1, U=1.0) = (r<1e-9) ? U : U*tanh(4*π*r/a)/(4*π*r/a)


# Functions with vector arguments
@scalar2vector CappedYukawa
@scalar2vector Hubbard
@scalar2vector Tanh

# Lattice operators
function gethubbard(lat, neighbors=[[0;0]]; mode=:nospin, format=:auto, kwargs...)
    """
    returns Dict(δL => Matrix(V(r_i-r_j+δL))_ij) where δL are vectors that
    connect unit cells. The set of δL's (in units of lattice vectors) is specified by 'neighbors'.
    """
    t(args...) = Hubbard(args...; kwargs...)
    ee_exchange = Hops(lat, neighbors, t; format=format)

    addspin(ee_exchange, mode)
end

function getshortrangedpotential(lat, V0, V1=0, V2=0; spin=true)

    dists, neighborhops = getneighborhops(lat;cellrange=4)
    hops1, hops2, hops3 = neighborhops[1:3]

    hops = V0 * hops1 

    if V1 != 0
        hops += V1 * hops2
    end
    if V2 != 0
        hops += V2 * hops3
    end

    if spin
        hops = kron(hops, ones(2,2))
        hops[zerokey(hops)][diagind(hops[zerokey(hops)])] .= 0
    end

    Hops(R=>complex(M) for (R,M) in hops)
end

# function getcappedyukawa(lat; cellrange=1, mode=:nospin, format=:auto, kwargs...)
#     t(args...) = CappedYukawa(args...; kwargs...)
#     ee_exchange = Hops(lat, neighbors, t; cellrange=cellrange, format=format)

#     addspin(ee_exchange, mode)
# end

function getcappedyukawa(lat, args...; spin=true, k0=1.0, a=nothing, U=1.0, kwargs...)
    t(args0...) = CappedYukawa(args0...; k0=k0, U=U)
    ee_exchange = Hops(lat, args..., t; kwargs...)

    if a!=nothing
        @info "Overwriting keyword argument k0 with a^(-1)."
        k0=1/a
    end

    if spin
        ee_exchange = kron(ee_exchange, ones(2,2))
        ee_exchange[zerokey(ee_exchange)][diagind(ee_exchange[zerokey(ee_exchange)])] .= 0
    end

    # addspin(ee_exchange, mode)
    ee_exchange
end

function gettanh(lat, args...; spin=true, a=1.0, U=1.0, kwargs...)
    t(args0...) = Tanh(args0...; a=a, U=U)
    ee_exchange = Hops(lat, args..., t; kwargs...)

    if spin
        ee_exchange = kron(ee_exchange, ones(2,2))
        ee_exchange[zerokey(ee_exchange)][diagind(ee_exchange[zerokey(ee_exchange)])] .= 0
    end

    # addspin(ee_exchange, mode)
    ee_exchange
end

# todo: implement https://en.wikipedia.org/wiki/Screened_Poisson_equation#Two_dimensions

# build_CappedYukawa(lat; mode=:nospin, format=:auto, kwargs...) = gethamiltonian(lat, r->CappedYukawa(r; kwargs...); mode=mode, format=format)
# gethamiltonianubbard(lat; mode=:nospin, format=:auto, kwargs...) = gethamiltonian(lat, r->Hubbard(r; kwargs...); mode=mode, format=format)


