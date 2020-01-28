using Distributed

"""
    General routines
    ================
    DOS = dos(...)
"""

function dos!(DOS, energy::Number, frequencies; broadening::Float64)
    DOS .+=(imag.( 1.0./(frequencies .- 1.0im .* broadening .- energy) ))
end
function dos!(DOS, energies::AbstractVector, args...; kwargs...)
    for ϵ in energies
        dos!(DOS, ϵ, args...; kwargs...)
    end
end

function dos_serial(h, ks::AbstractMatrix{Float64}, frequencies::AbstractVector{Float64}; Γ::Float64, kwargs...)
    L = size(ks)[2]
    DOS = zero(frequencies)

    ϵs = energies(h; kwargs...)

    @showprogress 1 "Computing DOS... " for k=eachcol(ks) # j=1:L
        dos!(DOS, ϵs(k), frequencies; broadening=Γ)
    end

    DOS / L / π
end

using SharedArrays

function dos_parallel(h, ks::AbstractMatrix{Float64}, frequencies::AbstractVector{Float64}; Γ::Float64, kwargs...)
    L = size(ks)[2]

    DOS = SharedArray(zero(frequencies))

    ϵs = energies(h; kwargs...)

    @sync @showprogress 1 "Computing DOS... " @distributed for j=1:L # over ks
        dos!(DOS, ϵs(ks[:,j]), frequencies; broadening=Γ)
    end

    DOS / L / π
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

using ..Utils: regulargrid

function dos(h, frequencies::AbstractVector{Float64}; klin::Int, kwargs...)
    ks = regulargrid(nk=klin^2)
    dos(h, ks, frequencies; kwargs...)
end

function dos(h, ks::AbstractMatrix{Float64}, frequencies::AbstractVector{Float64}; mode=:parallel, kwargs...)
    if mode==:parallel
        dos_parallel(h, ks, frequencies; kwargs...)
    elseif mode==:serial
        dos_serial(h, ks, frequencies; kwargs...)
    else
        error("Invalid mode value dos(...). Must be :parallel or :serial.")
    end
end

dos_dense(h, args...; kwargs...) = dos(h, args...; format=:dense, kwargs...)
dos_sparse(h, args...; kwargs...) = dos(h, args...; format=:sparse, kwargs...)


"""
    Alternative dos method with adaptive quadrature (better convergence)
    ====================================================================
    (DOS, ERROR) = dos_quad(...)
"""

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

using Cubature

function dos_quad_parallel(ϵs::Function, energies::AbstractVector{Float64};  Γ::Float64, kwargs...)
    fdim = length(energies)
    f0(k) = sum( imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) ) for ϵ=ϵs(k) )
    function f(ks, v)
        v[:] = pmap(f0, eachcol(ks))
    end
    xmin = [ 0.0, 0.0 ]
    xmax = [ 1.0, 1.0 ]

    (DOS, ERROR) = hcubature_v(fdim, f, xmin, xmax; kwargs...)

    DOS/π, ERROR
end