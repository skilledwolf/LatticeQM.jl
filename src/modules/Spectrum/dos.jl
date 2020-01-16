using Distributed

"""
    General routines
    ================
    DOS = dos(...)
"""

function dos_serial(ϵs::Function, ks::AbstractMatrix{Float64}, energies::AbstractVector{Float64}; Γ::Float64)
    L = size(ks)[2]
    DOS = zero(energies)

    for k=eachcol(ks) # j=1:L
        for ϵ in ϵs(k) #(ks[:,j])
            DOS .+= imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) )
        end
    end

    DOS / L / π
end

function dos_parallel(ϵs::Function, ks::AbstractMatrix{Float64}, energies::AbstractVector{Float64}; Γ::Float64)
    L = size(ks)[2]

    dos = @distributed (+) for j=1:L # over ks
        tmpdos = zero(energies)

        for ϵ in ϵs(ks[:,j]) # over bands at k
            tmpdos .+= imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) )
        end

        tmpdos / L / π
    end

    dos
end

function dos(ϵs::Function, ks::AbstractMatrix{Float64}, energies::AbstractVector{Float64}; Γ::Float64, mode=:parallel)
    if mode==:parallel
        dos_parallel(ϵs, ks, energies; Γ=Γ)
    elseif mode==:serial
        dos_serial(ϵs, ks, energies; Γ=Γ)
    else
        error("Invalid mode value dos(...). Must be :parallel or :serial.")
    end
end


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


"""
    Interface to dense/sparse routines
    ==================================
"""

function dos_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, ϵ::AbstractVector{Float64}; Γ::Float64, mode=:parallel)
    dos(energies(hamiltonian, format=:dense), ks, ϵ; Γ=Γ, mode=mode)
end

function dos_sparse(hamiltonian::Function, ks::AbstractMatrix{Float64}, ϵ::AbstractVector{Float64}; Γ::Float64, mode=:parallel, kwargs...)
    dos(energies(hamiltonian; format=:sparse, kwargs...), ks, ϵ; Γ=Γ, mode=mode)
end

# function dos_dense_parallel(hamiltonian::Function, ks::AbstractMatrix{Float64}, energies::AbstractVector{Float64}; Γ::Float64)
#
#     L = size(ks)[2]
#
#     dos = @distributed (+) for j=1:L # over ks
#         tmpdos = zero(energies)
#
#         for ϵ in ϵs_dense(hamiltonian)(ks[:,j]) # over bands at k
#             tmpdos .+= imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) )
#         end
#
#         tmpdos / L / π
#     end
#
#     dos
# end

# function dos_dense_parallel(hamiltonian::Function, ks::AbstractMatrix{Float64}, ωmin::Float64, ωmax::Float64, num::Int)
#
#     energies = collect(range(ωmin, length=num, stop=ωmax))
#     Γ = 0.1 * abs(ωmax-ωmin)/(num-1)
#
#     energies, dos_dense_parallel(hamiltonian, ks, energies; Γ=Γ)
# end
