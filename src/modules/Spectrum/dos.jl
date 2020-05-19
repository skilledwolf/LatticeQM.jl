using Distributed
using SharedArrays
using ..Utils: regulargrid


"""
    getdos(h, ωs; klin, Γ, kwargs...)

Computes the density of states of operator h(k) on the entire Brillouine zone,
discretized on a grid with \$ k_{lin} \\times k_{lin} \$ points. 
and for the frequencies ωs=(ω1, ω2, ...). The paremter \$\\Gamma\$ is the energy broadening.

Accepts the same kwargs as dos(h, ks, ω).

Note: the current implementation only works for a two-dimensional Brillouine zone.
Might change in the future, but for now use dos(h, ks, ω; Γ) syntax if needed.
"""
getdos(h, ωs; kwargs...) = (DOS=zero(ωs); getdos!(DOS, h, ωs; kwargs...))
function getdos!(DOS, h, ωs::AbstractVector; klin::Int, kwargs...)
    ks = regulargrid(nk=klin^2)
    getdos!(DOS, h, ks, ωs; kwargs...)
end


"""
    getdos(h, ks, ω; Γ, mode=:parallel, format=:auto)

Computes the density of states of operator h(k) using the points ks=(k1,k2,...)
and for the frequencies ω=(ω1, ω2, ...). The paremter \$\\Gamma\$ is the energy broadening.

Mode can be :parallel or :serial, format can be :auto, :sparse or :dense.
"""
getdos(h, ks, ωs; kwargs...) = (DOS=zero(ωs); dos!(DOS, h, ks, ωs; kwargs...))
function getdos!(DOS, h, ks::AbstractMatrix{<:Real}, frequencies::AbstractVector{<:Number}; parallel=true, kwargs...)
    if (parallel && nprocs()<2)
        parallel=false
    end
    
    parallel ? dos_parallel!(DOS, h, ks, frequencies; kwargs...) : dos_serial!(DOS, h, ks, frequencies; kwargs...)
end

getdos_dense(h, args...; kwargs...) = getdos(h, args...; format=:dense, kwargs...)
getdos_sparse(h, args...; kwargs...) = getdos(h, args...; format=:sparse, kwargs...)


dos!(DOS, energies::AbstractVector, ωs; kwargs...) = foreach(x->dos!(DOS, x, ωs; kwargs...), energies)
function dos!(DOS, energy::Number, ωs; broadening::Number)
    DOS .+= imag.( 1.0./(ωs .- 1.0im .* broadening .- energy) )
end

function dos_serial!(DOS, h, ks::AbstractMatrix{<:Real}, frequencies::AbstractVector{<:Number}; Γ::Number, kwargs...)
    L = size(ks)[2]
    ϵs = energies(h; kwargs...)

    @showprogress 6 "Computing DOS... " for k=eachcol(ks) # j=1:L
        dos!(DOS, ϵs(k), frequencies; broadening=Γ)
    end
    DOS ./= (L / π)

    DOS
end

function dos_parallel!(DOS, h, ks::AbstractMatrix{<:Real}, frequencies::AbstractVector{<:Number}; Γ::Number, kwargs...)
    L = size(ks)[2]
    ϵs = energies(h; kwargs...)

    DOS0 = SharedArray(Array(DOS))
    @sync @showprogress 6 "Computing DOS... " @distributed for j=1:L # over ks
        dos!(DOS0, ϵs(ks[:,j]), frequencies; broadening=Γ)
    end
    DOS[:] .= DOS0[:]
    DOS[:] ./= (L / π)

    DOS
end



# function dos_parallel(h, ks::AbstractMatrix{Float64}, frequencies::AbstractVector{Float64}; Γ::Float64, kwargs...)
#     L = size(ks)[2]
#
#     ϵs = energies(h; kwargs...)
#
#     dos = @sync @showprogress 1 "Computing DOS... " @distributed (+) for j=1:L # over ks
#         tmpdos = zero(frequencies)
#         dos!(tmpdos, ϵs(ks[:,j]), frequencies; broadening=Γ)
#         tmpdos
#     end
#
#     dos / L / π
# end

# using HCubature
#
# function dos_quad(ϵs::Function, energies::AbstractVector{Float64};  Γ::Float64, kwargs...)
#     fdim = length(energies)
#     f(k) = sum( imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) ) for ϵ=ϵs(k) )
#     xmin = [ 0.0, 0.0 ]
#     xmax = [ 1.0, 1.0 ]
#
#     (DOS, ERROR) = hcubature(f, xmin, xmax; kwargs...)
#
#     DOS/π, ERROR
# end

# using Cubature

# function dos_quad_parallel(ϵs::Function, energies::AbstractVector{T};  Γ::T, kwargs...) where {T<:Number}
#     fdim = length(energies)
#     f0(k) = sum( imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) ) for ϵ=ϵs(k) )
#     function f(ks, v)
#         v[:] = pmap(f0, eachcol(ks))
#     end
#     xmin = [ 0.0, 0.0 ]
#     xmax = [ 1.0, 1.0 ]

#     (DOS, ERROR) = hcubature_v(fdim, f, xmin, xmax; kwargs...)

#     DOS/π, ERROR
# end