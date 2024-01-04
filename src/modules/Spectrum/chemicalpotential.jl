
##########################################################################################
# Low-level chemicalpotential solver
##########################################################################################

chemicalpotential(bands::AbstractArray, kgrid, filling::Real; T::Real = 0.0, kwargs...) = (T == zero(T)) ? chemicalpotential_0(bands, filling; kwargs...) : chemicalpotential_T(bands, kgrid, filling; T = T, kwargs...)

chemicalpotential_0(energies, args...; kwargs...) = (en = energies[:]; chemicalpotential_0!(en, args...; kwargs...))
function chemicalpotential_0!(energies::AbstractVector{<:Real}, filling::Real)

    if filling == 1.0
        en = maximum(energies)
    else
        i = 1 + floor(Int, filling * (length(energies) - 1)) # fill a fraction of states according to ,,filling''
        en = sum(partialsort!(energies, i:i+1)) / 2 # oldversion: sort!(energies); e1, e2 = energies[i:i+1]
    end

    en
end


# import NLsolve # add 1-2 seconds of load time to the package :(
import LatticeQM.Utils
import LatticeQM.Utils: fermidirac

import Roots

function chemicalpotential_T(bands, kgrid, filling::T1; T::T1 = 0.01) where {T1<:Real}

    #  solve  n(μ*) = n_occ for μ*
    # μ0 = chemicalpotential_0(bands, filling) # initial guess
    # μ, res = brentq(x -> getfilling(bands, kgrid; μ = x) - filling, μ0 - 2 * T, μ0 + 10 * T; full_output = true)
    μ = Roots.find_zero(x -> getfilling(bands, kgrid; μ = x) - filling, (minimum(bands), maximum(bands)))

    # if !res.converged
    #     println("WARNING: Calculation of chemical potential for finite T did not converge.")
    #     println("         Using chemical potential for T=0 instead (=Fermi energy).")
    #     μ = chemicalpotential_0(bands, filling)
    # end

    μ
end


##########################################################################################
# Interface of chemicalpotential to Hamiltonian type
##########################################################################################

function chemicalpotential(H, ks, filling::Real; multimode = :distributed, kwargs...)
    H = Utils.getelectronsector(H)
    chemicalpotential(bandmatrix(H, ks; multimode = multimode, progress_label="Chemical potential")[1], ks, filling; kwargs...)
end

function chemicalpotential(H, ks, fillings::AbstractVector; multimode = :distributed, kwargs...)
    H = Utils.getelectronsector(H)
    bands = bandmatrix(H, ks; multimode = multimode)[1] # compute once to determine multiple fillings later on
    [chemicalpotential(bands, ks, filling; kwargs...) for filling in fillings]
end


import ..Structure: Mesh, meshweights


##########################################################################################
# getfilling for non-regular grids
##########################################################################################
import Statistics

getfillings(bands::AbstractMatrix, μs::AbstractVector, args...; kwargs...) = [getfilling(bands, args...; μ = μ, kwargs...) for μ in μs]

getfilling(bands::AbstractMatrix, kmesh::Mesh; kwargs...) = getfilling(bands, meshweights(kmesh); kwargs...)

getfilling(bands::AbstractArray, kpoints::AbstractMatrix; kwargs...) = Statistics.mean(fermidirac(bands; kwargs...))
getfilling(bands::AbstractArray; kwargs...) = Statistics.mean(fermidirac(bands; kwargs...))
getfilling(bands::AbstractMatrix, weights::AbstractVector; kwargs...) = sum(weights' .* fermidirac(bands; kwargs...)) / size(bands, 1)



##########################################################################################
# Interface of getfilling to Hamiltonian type
##########################################################################################

function filling(H, ks, μ; multimode = :distributed, kwargs...)

    getfilling(bandmatrix(Utils.getelectronsector(H), ks; multimode = multimode)[1], ks; μ = μ, kwargs...)
end

import ..Structure: regulargrid

function filling(H, μ; multimode = :distributed, nk = 10, kwargs...)
    ks = regulargrid(nk = nk^2)
    getfilling(bandmatrix(Utils.getelectronsector(H), ks; multimode = multimode)[1]; μ = μ, kwargs...)
end

function fillings(H, ks, μs; kwargs...)
    bandmatrix = bandmatrix(Utils.getelectronsector(H), ks; multimode = multimode)[1]

    getfillings(bandmatrix, μs, ks; kwargs...)
end

function fillings(H, μs; nk=10, kwargs...)
    ks = regulargrid(nk = nk^2)
    bandmatrix = bandmatrix(Utils.getelectronsector(H), ks; multimode = multimode)[1]

    getfillings(bandmatrix, μs, ks; kwargs...)
end



##########################################################################################
# Alternate all-julia implementation: (requires NLsolve, which loads slowly)
##########################################################################################

# function chemicalpotential_T(bands, kgrid, filling::T1; T::T1 = 0.01) where {T1<:Real}

#     # Initial guess for T=0
#     μ0 = chemicalpotential_0(bands, filling)

#     ##########################################
#     ###  Solve  n(μ*) = N_occ for μ*
#     ##########################################

#     n(μ::AbstractFloat) = getfilling(bands, kgrid; μ = μ)

#     function δn!(δn::T, μ::T) where {T2<:AbstractFloat, T<:AbstractArray{T2}}
#         δn[1] = n(μ[1]) - filling
#         nothing
#     end
#     sol = NLsolve.nlsolve(δn!, [μ0])#; ftol=1e-6, iterations=1500)
#     if !NLsolve.converged(sol)
#         println("WARNING: Calculation of chemical potential for finite T did not converge.")
#         println("         Using chemical potential for T=0 instead (=Fermi energy).")
#         μ = μ0
#     else
#         μ = sol.zero[1]
#     end

#     μ
# end
