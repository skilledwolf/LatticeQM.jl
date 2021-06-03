# TODO: implement chemicalpotential!(bands::AbstractMatrix, ...)
# the repeated allocation of bands is problematic for large systems!

getelectronsector(H::Function) = H

function chemicalpotential(H, ks, filling::Real; multimode=:distributed, kwargs...)

    chemicalpotential(bandmatrix(getelectronsector(H), ks; multimode=multimode), filling; kwargs...)
end

function chemicalpotential(bands::AbstractMatrix, filling::Real; kwargs...)
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

function chemicalpotential!(energies::AbstractVector{<:Real}, nk::Int, filling::Float64; T::Real=0.0, kwargs...)
    if T==zero(T)
        return chemicalpotential_0!(energies, filling; kwargs...)
    else
        return chemicalpotential_T!(energies, nk, filling; T=T, kwargs...)
    end
end

function chemicalpotential_0!(energies::AbstractVector{<:Real}, filling::Real)

    if filling == 1.0
        return maximum(energies)
    end

    i = 1 + floor(Int, filling * (length(energies)-1)) # fill a fraction of states according to ,,filling''

    e1, e2 = partialsort!(energies, i:i+1) # oldversion: sort!(energies); e1, e2 = energies[i:i+1]

    (e1+e2)/2
end


# import NLsolve # add 1-2 seconds of load time to the package :(
import ..Utils: fermidirac, brentq

function chemicalpotential_T!(energies::AbstractVector{T1}, nk::Int, filling::T1; T::T1=0.01) where T1<:Real

    d = div(length(energies),nk) # number of bands

    # Initial guess for T=0
    μ0 = chemicalpotential_0!(energies, filling)

    ##########################################
    ###  Solve  n(μ*) = N_occ for μ*
    ##########################################

    n(μ::AbstractFloat) = sum(fermidirac(ϵ-μ; T=T) for ϵ=energies)/nk

    μ, res = brentq(x->n(x)-d*filling, μ0-2*T, μ0+10*T; full_output=true)
    if !res.converged
        println("WARNING: Calculation of chemical potential for finite T did not converge.")
        println("         Using chemical potential for T=0 instead (=Fermi energy).")
        μ = μ0
    end

    # function δn!(δn::T, μ::T) where {T2<:AbstractFloat, T<:AbstractArray{T2}}
    #     δn[1] = n(μ[1]) - d * filling
    #     nothing
    # end
    # sol = NLsolve.nlsolve(δn!, [μ0])#; ftol=1e-6, iterations=1500)
    # if !NLsolve.converged(sol)
    #     println("WARNING: Calculation of chemical potential for finite T did not converge.")
    #     println("         Using chemical potential for T=0 instead (=Fermi energy).")
    #     μ = μ0
    # else
    #     μ = sol.zero[1]
    # end
    
    μ
end
