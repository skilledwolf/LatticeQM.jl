function setfilling!(hops, lat, filling; nk=100, kwargs...)
    kgrid = regulargrid(nk=nk)
    setfilling!(hops, lat, kgrid, filling; kwargs...)
end

function setfilling!(hops, lat, kgrid, filling; T=0.0)
    hops0 = DenseHops(hops)
    μ = chemicalpotential(getbloch(hops0), kgrid, filling; T=T)
    addchemicalpotential!(hops, lat, -μ)

    μ
end

function addchemicalpotential!(hops, lat::Lattice, μ::T; localdim::Int=-1) where T<:AbstractVector{<:Float64}
    zero0 = zeros(Int, latticedim(lat))
    N = countorbitals(lat)

    if localdim < 0 # if localdim is not set, we determine it from matrix dimensions
        D = hopdim(hops)
        d = div(D, N)
    else
        d = localdim
    end

    @assert N == length(μ)

    newhops = Dict( zero0 => sparse( (1.0+0.0im).* Diagonal(kron(μ, ones(d))) ) )
    addhops!(hops, newhops)

    nothing
end

function addchemicalpotential!(hops, lat::Lattice, μ::Function; kwargs...)
    R = allpositions(lat)
    addchemicalpotential!(hops, lat, [μ(r) for r=eachcol(R)]; kwargs...)

    nothing
end
addchemicalpotential!(hops, lat::Lattice, μ::Float64; kwargs...) = addchemicalpotential!(hops, lat, μ.*ones(countorbitals(lat)); kwargs...)

@legacyalias addinterlayerbias! add_interlayerbias!
@legacyalias addinterlayerbias! add_transversepotential!
@legacyalias addinterlayerbias! addtransversepotential!
function addinterlayerbias!(hops, lat::Lattice, V::Float64; d=3.0, kwargs...)

    # Only go through the trouble of constructing this matrix for finite V
    if abs(V) ≈ 0
        return nothing
    end

    # Get z coordinates and scale them into the unit range
    layer = extrapositions(lat, "z")
    min = minimum(layer); max = maximum(layer)

    if abs(min-max) ≈ 0
        error("Requires at least two layers.")
    end

    layer .= (layer .- min)./(max-min)
    μ = V .* (layer .- 0.5)

    add_chemicalpotential!(hops, lat, μ)

    nothing
end

###################################################################################################
# Backwards compatibility
###################################################################################################
export set_filling!
@legacyalias setfilling! set_filling!

export add_chemicalpotential!
@legacyalias addchemicalpotential! add_chemicalpotential!
