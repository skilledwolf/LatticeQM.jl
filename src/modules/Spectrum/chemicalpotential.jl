
##########################################################################################
# Low-level chemicalpotential solver
##########################################################################################

import LatticeQM.Utils
import LatticeQM.Utils: fermidirac
import ..Structure: Mesh, meshweights
import Roots

function chemicalpotential(bands::AbstractArray, kgrid, filling::Real; T::Real=0.0)
    0 <= filling <= 1 || throw(ArgumentError("filling must lie in [0, 1]"))
    T >= 0 || throw(ArgumentError("temperature must be nonnegative"))
    T == 0 ? chemicalpotential_0(bands, kgrid, filling) :
             chemicalpotential_T(bands, kgrid, filling; T)
end

chemicalpotential_0(bands, kgrid, filling::Real) = chemicalpotential_0(bands, filling)
chemicalpotential_0(bands, mesh::Mesh, filling::Real) = chemicalpotential_0(bands, meshweights(mesh), filling)
chemicalpotential_0(energies, filling::Real) = chemicalpotential_0!(vec(copy(energies)), filling)

# At T=0 a finite grid cannot represent arbitrary fillings. Return a gap
# midpoint when the target is attainable, otherwise the energy of the
# partially occupied level (fermidirac assigns half occupation there).
function _chemicalpotential_sorted(energies, weights, filling)
    0 <= filling <= 1 || throw(ArgumentError("filling must lie in [0, 1]"))
    filling == 0 && return prevfloat(float(first(energies)))
    filling == 1 && return nextfloat(float(last(energies)))
    occupied, correction = 0.0, 0.0
    for i in eachindex(energies)
        # Compensated accumulation keeps a large uniform mesh from missing
        # an exactly attainable filling through summation roundoff.
        increment = weights[i] - correction
        total = occupied + increment
        correction = (total - occupied) - increment
        occupied = total
        # All copies of a degenerate level must be considered together.
        i < length(energies) && energies[i] == energies[i+1] && continue
        if isapprox(occupied, filling; rtol=8eps(Float64), atol=8eps(Float64))
            return i < length(energies) ? energies[i]/2 + energies[i+1]/2 : energies[i]
        elseif occupied > filling
            return energies[i]
        end
    end
    last(energies)
end

function chemicalpotential_0!(energies::AbstractVector{<:Real}, filling::Real)
    0 <= filling <= 1 || throw(ArgumentError("filling must lie in [0, 1]"))
    filling == 0 && return prevfloat(float(minimum(energies)))
    filling == 1 && return nextfloat(float(maximum(energies)))
    occupied = filling * length(energies)
    count = round(Int, occupied)
    if 0 < count < length(energies) && isapprox(occupied, count; rtol=8eps(Float64))
        edges = partialsort!(energies, count:count+1)
        return edges[1]/2 + edges[2]/2
    end
    partialsort!(energies, ceil(Int, occupied))
end

function chemicalpotential_0(bands::AbstractMatrix, weights::AbstractVector, filling::Real)
    length(weights) == size(bands, 2) || throw(DimensionMismatch("one weight is required per k point"))
    all(>=(0), weights) && isapprox(sum(weights), 1) ||
        throw(ArgumentError("k-point weights must be nonnegative and sum to one"))
    energies = vec(bands)
    order = sortperm(energies)
    stateweights = repeat(weights ./ size(bands, 1); inner=size(bands, 1))
    _chemicalpotential_sorted(energies[order], stateweights[order], filling)
end

function chemicalpotential_T(bands, kgrid, filling::Real; T::Real=0.01)
    filling == 0 && return -Inf
    filling == 1 && return Inf
    # At μ=Emin+T*logit(f), every state's occupation is at most f;
    # the analogous Emax bound gives at least f, including extreme fillings.
    shift = T * (log(filling) - log1p(-filling))
    lo, hi = extrema(bands)
    lo == hi && return lo + shift
    Roots.find_zero(μ -> getfilling(bands, kgrid; μ, T) - filling,
                    (lo + shift - T, hi + shift + T))
end


##########################################################################################
# Interface of chemicalpotential to Hamiltonian type
##########################################################################################

function chemicalpotential(H, ks, filling::Real; multimode = :distributed, hidebar=false, kwargs...)
    H = Utils.getelectronsector(H)
    chemicalpotential(bandmatrix(H, ks; multimode=multimode, hidebar=hidebar,
                                  progress_label="Chemical potential")[1],
                      ks, filling; kwargs...)
end

function chemicalpotential(H, ks, fillings::AbstractVector; multimode = :distributed, hidebar=false, kwargs...)
    H = Utils.getelectronsector(H)
    bands = bandmatrix(H, ks; multimode=multimode, hidebar=hidebar,
                       progress_label="Chemical potential")[1]
    [chemicalpotential(bands, ks, filling; kwargs...) for filling in fillings]
end

##########################################################################################
# getfilling for non-regular grids
##########################################################################################
import Statistics

getfilling(bands::AbstractMatrix, kmesh::Mesh; kwargs...) = getfilling(bands, meshweights(kmesh); kwargs...)

getfilling(bands::AbstractArray, kpoints::Union{AbstractMatrix,DiscretePath}; kwargs...) = Statistics.mean(fermidirac(bands; kwargs...))
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

