function setfilling!(H, lat, filling; nk=100, kwargs...)
    kgrid = regulargrid(nk=nk)
    setfilling!(H, lat, kgrid, filling; kwargs...)
end

function setfilling!(H, lat, kgrid, filling; T=0.0)
    μ = chemicalpotential(H, kgrid, filling; T=T)
    addchemicalpotential!(H, lat, -μ)
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
    layer = positions(lat, 3)
    min = minimum(layer); max = maximum(layer)

    if isapprox(abs(min-max), 0; atol=sqrt(eps()))
        error("Requires at least two layers.")
    end

    layer .= (layer .- min)./(max-min)
    μ = V.(eachcol(R)) .* (layer .- 0.5)

    addchemicalpotential!(hops, lat, vec(μ))

    nothing
end