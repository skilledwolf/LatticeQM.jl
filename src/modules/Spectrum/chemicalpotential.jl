# TODO: implement chemicalpotential!(bands::AbstractMatrix, ...)
# the repeated allocation of bands is problematic for large systems!

using ..TightBinding: getelectronsector

function chemicalpotential(H, ks, filling::Float64; multimode=:distributed, kwargs...)

    chemicalpotential(bandmatrix(getelectronsector(H), ks; multimode=multimode), filling; kwargs...)
end

function chemicalpotential(bands::AbstractMatrix, filling::Float64; kwargs...)
    en = bands[:]
    nk = size(bands,2)

    chemicalpotential!(en, nk, filling; kwargs...)
end

# function chemicalpotential(energies::AbstractVector, ks::AbstractMatrix, filling::Float64; kwargs...)
#     en = energies[:]
#     nk = size(ks,2)
#
#     chemicalpotential!(en, nk, filling; kwargs...)
# end

function chemicalpotential!(energies::AbstractVector{T1}, nk::Int, filling::Float64; T::T1=0.0) where T1<:Real
    if T==zero(T)
        return chemicalpotential_0!(energies, filling)
    else
        return chemicalpotential_T!(energies, nk, filling; T=T)
    end
end

function chemicalpotential_0!(energies::AbstractVector{T1}, filling::T1) where T1<:Real

    i = floor(Int, filling * length(energies)) # fill a fraction of states according to ,,filling''
    e1, e2 = partialsort!(energies, i:i+1) # oldversion: sort!(energies); e1, e2 = energies[i:i+1]

    (e1+e2)/2
end


import NLsolve: nlsolve, converged # add 2-3 seconds of load time to the package :(
import ..Utils: fermidirac

function chemicalpotential_T!(energies::AbstractVector{T1}, nk::Int, filling::T1; T::T1=0.01) where T1<:Real

    d = div(length(energies),nk) # number of bands

    # Initial guess for T=0
    μ0 = chemicalpotential_0!(energies, filling)

    ##########################################
    ###  Solve  n(μ*) = N_occ for μ*
    ##########################################

    n(μ::AbstractFloat) = sum(fermidirac(ϵ-μ; T=T) for ϵ=energies)/nk

    function δn!(δn::T, μ::T) where {T2<:AbstractFloat, T<:AbstractArray{T2}}
        δn[1] = n(μ[1]) - d * filling
        nothing
    end

    sol = nlsolve(δn!, [μ0])#; ftol=1e-6, iterations=1500)
    # @assert converged(sol) "Could not calculate the chemical potential for finite T."
    if !converged(sol)
        println("WARNING: Calculation of chemical potential for finite T did not converge.")
        println("         Using chemical potential for T=0 instead (=Fermi energy).")
        μ = μ0
    else
        μ = sol.zero[1]
    end
    
    μ
end
