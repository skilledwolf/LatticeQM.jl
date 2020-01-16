# TODO: implement chemicalpotential!(bands::AbstractMatrix, ...)
# the repeated allocation of bands is problematic for large systems!

chemicalpotential(hops::AnyHops, args...; kwargs...) = chemicalpotential(getbloch(hops), args...; kwargs...)

function chemicalpotential(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64; kwargs...)

    chemicalpotential(bandmatrix(hamiltonian, ks)[:], ks, filling; kwargs...)
end

function chemicalpotential(energies::AbstractVector, ks::AbstractMatrix, filling::Float64; kwargs...)

    en = energies[:]
    nk = size(ks,2)

    chemicalpotential!(en, nk, filling; kwargs...)
end

function chemicalpotential!(energies::AbstractVector{Float64}, nk::Int, filling::Float64; T::AbstractFloat=0.0)
    if T==0.0
        return chemicalpotential_0!(energies, filling)
    else
        return chemicalpotential_T!(energies, nk, filling; T=T)
    end
end

function chemicalpotential_0!(energies::AbstractVector{Float64}, filling::Float64)

    i = floor(Int, filling * length(energies)) # fill a fraction of states according to ,,filling''
    e1, e2 = partialsort!(energies, i:i+1) # oldversion: sort!(energies); e1, e2 = energies[i:i+1]

    (e1+e2)/2
end

using NLsolve

using ..Utils: fermidirac

function chemicalpotential_T!(energies::AbstractVector{Float64}, nk::Int, filling::Float64; T::Float64=0.01)

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

    sol = nlsolve(δn!, [μ0])
    @assert converged(sol)
    μ = sol.zero[1]

    μ
end

###################################################################################################
# Backwards compatibility
###################################################################################################
@legacyalias chemicalpotential chemical_potential

