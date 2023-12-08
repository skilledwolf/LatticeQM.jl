import ..Structure: regulargrid
import ..Structure.Lattices: Lattice, latticedim, countorbitals, allpositions, positions
import ..TightBinding: Hops, AbstractHops, zerokey, hopdim, addhops!

import LinearAlgebra: Diagonal

import SparseArrays: sparse

function setfilling!(H, filling; nk=100, kwargs...)
    kgrid = regulargrid(nk=nk)
    setfilling!(H, kgrid, filling; kwargs...)
end

import ..Utils
import ..Spectrum

function setfilling!(H, kgrid, filling; kwargs...)
    μ = Spectrum.chemicalpotential(Utils.dense(H), kgrid, filling; kwargs...)
    addchemicalpotential!(H, -μ)
    μ
end

function addchemicalpotential!(hops::Hops, μ::Real)
    inds = diagind(hops[zerokey(hops)])
    hops[zerokey(hops)][inds] .+= μ
    hops
end

function addchemicalpotential!(hops::Hops, lat::Lattice, μ::T; localdim::Int=-1) where T<:AbstractVector{<:Float64}
    # zero0 = zeros(Int, latticedim(lat))
    N = countorbitals(lat)

    if localdim < 0 # if localdim is not set, we determine it from matrix dimensions
        D = hopdim(hops)
        d = div(D, N)
    else
        d = localdim
    end

    @assert N == length(μ)

    newhops = Hops( zerokey(hops) => sparse( (1.0+0.0im).* Diagonal(kron(μ, ones(d))) ) )
    addhops!(hops, newhops)

    nothing
end

function addchemicalpotential!(hops, lat::Lattice, μ::Function; kwargs...)
    R = allpositions(lat)
    addchemicalpotential!(hops, lat, [μ(r) for r=eachcol(R)]; kwargs...)

    nothing
end
addchemicalpotential!(hops, lat::Lattice, μ::Float64; kwargs...) = addchemicalpotential!(hops, lat, μ.*ones(countorbitals(lat)); kwargs...)

function addinterlayerbias!(hops, lat, V::Real; kwargs...)
    # Only go through the trouble of constructing this matrix for finite V
    if isapprox(V,0; atol=sqrt(eps()))
        return nothing
    end
    addinterlayerbias!(hops,lat,x->V; kwargs...)
end
function addinterlayerbias!(hops, lat::Lattice, V::Function; d=3.0, kwargs...)

    # Get z coordinates and scale them into the unit range
    R = allpositions(lat)
    Z = positions(lat, 3)
    min = minimum(Z); max = maximum(Z)

    if isapprox(abs(min-max), 0; atol=sqrt(eps()))
        error("Requires at least two layers.")
    end

    Z .= (Z .- min)./(max-min)
    μ = V.(eachcol(R)) .* (Z .- 0.5)

    addchemicalpotential!(hops, lat, vec(μ); kwargs...)

    nothing
end