using SparseArrays

function unit(i::Int,j::Int, N::Int)
    mat = zeros(ComplexF64, N, N)
    mat[i,j] = 1.0

    mat
end

####################################################################################################
#  Chemical potential
####################################################################################################

function chemicalpotential!(hops, lat::Lattice, μ::T; localdim::Int=1) where T<:AbstractVector{<:Float64}

    zero0 = zeros(Int, lattice_dim(lat))

    N = atom_count(lat)
    @assert N == length(μ)

    newhops = Dict( zero0 => sparse( (1.0+0.0im).* Diagonal(kron(μ, ones(localdim))) ) )
    add_hoppings!(hops, newhops)

    nothing
end

function chemicalpotential!(hops, lat::Lattice, μ::Function; kwargs...)
    R = positionsND(lat)
    chemicalpotential!(hops, lat, [μ(r) for r=eachcol(R)]; kwargs...)

    nothing
end
chemicalpotential!(hops, lat::Lattice, μ::Float64; kwargs...) = chemicalpotential!(hops, lat, μ.*ones(atom_count(lat)); kwargs...)

function transversepotential!(hops, lat::Lattice, V::Float64; d=3.0, kwargs...)

    # Only go through the trouble of constructing this matrix for finite V
    if abs(V) ≈ 0
        return nothing
    end

    # Get z coordinates and scale them into the unit range
    layer = get_positions_in(lat, "z")
    min = minimum(layer); max = maximum(layer)

    if abs(min-max) ≈ 0
        error("Requires at least two layers.")
    end

    layer .= (layer .- min)./(max-min)
    μ = V .* (layer .- 0.5)

    chemicalpotential!(hops, lat, μ)

    nothing
end

####################################################################################################
#  Zeeman fields
####################################################################################################

function zeeman_hops(args...; kwargs...)
    newhops = Dict{Vector{Int},SparseMatrixCSC{ComplexF64}}()
    zeeman!(newhops, args...; kwargs...)

    newhops
end

function zeeman!(hops, lat::Lattice, Mv::Function)
    zero0 = zeros(Int, lattice_dim(lat))

    N = atom_count(lat)
    R = positionsND(lat)

    σn(vec) = sum(vec[i] .* σs[i] for i=1:3)

    newhops = Dict(zero0 => sum(kron(unit(i,i,N), σn(Mv(R[:,i]))) for i=1:N))
    add_hoppings!(hops, newhops)

    nothing
end

function zeeman!(hops, lat::Lattice, Mv::Vector{Float64}; format=:dense)
    # Only go through the trouble of constructing this matrix for finite Mv
    if norm(Mv) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))

    N = atom_count(lat)
    R = positions(lat)

    σn = sum(Mv[i] .* σs[i] for i=1:3)

    newhops = Dict( zero0 => kron(Matrix(1.0I, N, N), σn) )
    add_hoppings!(hops, newhops)

    nothing
end
zeeman!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)

function zeeman_staggered!(hops, lat::Lattice, Mv::Vector{Float64}; format=:auto)
    """
    Generate a Zeeman term in the Hamiltonian given a Magnetization  M(r),
    where r is a position vector.

    lat::Lattice      Lattice container
    M::Function       returns magnetization M(r)::Vector{Float} at site r
    """

    # Only go through the trouble of constructing this matrix for finite Mv
    if norm(Mv) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "sublattice")

    N = atom_count(lat)
    R = positions(lat)
    sublattice = (get_positions_in(lat, "sublattice") .- 0.5) # +Mz/2 on A and -Mz/2 on B

    σn = sum(Mv[i] .* σs[i] for i=1:3)

    hops[zero0] += sum(kron(unit(i,i,N), sublattice[i] .* σn) for i=1:N)

    nothing
end
zeeman_staggered!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_staggered!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)


function zeeman_staggered_noncol!(hops, lat::Lattice, M0::Float64, ϕ=0.0::Float64; format=:auto)
    """
    Generate a Zeeman term in the Hamiltonian given a Magnetization  M(r),
    where r is a position vector.

    lat::Lattice      Lattice container
    M::Function       returns magnetization M(r)::Vector{Float} at site r
    """

    # Only go through the trouble of constructing this matrix for finite Mv
    if abs(M0) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "sublattice")

    N = atom_count(lat)
    R = positions(lat)
    sublattice = get_positions_in(lat, "sublattice") # +Mz/2 on A and -Mz/2 on B

    σn(Mv::T) where {T<:AbstractVector{Float64}} = sum(Mv[i] .* σs[i] for i=1:3)
    Mv(index) = [0.0,0.0,1.0] * index + [cos(ϕ), sin(ϕ), .0] * (1.0-index)

    hops[zero0] += sum(kron(unit(i,i,N), M0 .* σn(Mv(sublattice[i]))) for i=1:N)

    nothing
end


function zeeman_layered_noncol!(hops, lat::Lattice, M1::Vector{Float64}, M2::Vector{Float64})
    """
    Generate a Zeeman term in the Hamiltonian given a with a Magnetization vector
    that is
    """

    # Only go through the trouble of constructing this matrix for finite Mv
    if norm(M1) ≈ 0 && norm(M2) ≈ 0
        return 0.0
    end

    zero0 = zeros(Int, lattice_dim(lat))
    assert_dimension(lat, "z")

    N = atom_count(lat)
    R = positions(lat)

    # Get z coordinates and scale them into the unit range
    layer = get_positions_in(lat, "z")
    min = minimum(layer)
    max = maximum(layer)

    if abs(min-max) ≈ 0
        error("Requires at least two layers.")
    end

    layer .= (layer .- min)./(max-min)

    # Define the magnetization as function of direction
    Mv(z) = M1 .* (1.0-z) + M2 .* (z)
    σn(z) = sum(Mv(z)[i] .* σs[i] for i=1:3)

    # hops[zero0] += blockdiag([sparse(σn(z_scaled)) for z_scaled=layer]...)
    hops[zero0] += sum(kron(unit(i,i,N), σn(layer[i])) for i=1:N)

    nothing
end
zeeman_layered_noncol!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_layered_noncol!(hops, lat::Lattice, [M0,0.0,0.0], [0.0,0.0,M0]; kwargs...)
zeeman_layered!(hops, lat::Lattice, M0::Float64; kwargs...) = zeeman_layered_noncol!(hops, lat::Lattice, [0.0,0.0,M0], [0.0,0.0,-M0]; kwargs...)#zeeman_field_layered!(hops, lat::Lattice, [0.0,0.0,M0]; kwargs...)

