
mutable struct Mesh{A<:AbstractMatrix,B<:AbstractVector,C<:AbstractVector{<:Int}}
    points::A # each column is a point
    measure::B # each point has a measure of it's corresponding plaquette
    multiplicity::C # in case we have to keep track of a hidden multiplicity
end

function Mesh(points::AbstractMatrix, measure::AbstractVector)
    multiplicity = ones(Int, size(points,2))
    Mesh(points, measure, multiplicity)
end

function Mesh(points::AbstractMatrix)
    L = size(points,2)
    measure = 1/L * ones(size(points,2))
    multiplicity = ones(Int, size(points, 2))
    Mesh(points, measure, multiplicity)
end

meshweights(data::Mesh) = data.measure .* data.multiplicity

Base.length(data::Mesh) = Base.size(data.points, 2) # columns in points
Base.iterate(data::Mesh) = length(data) > 0 ? (data[Base.firstindex(data)], Base.firstindex(data)) : nothing
Base.iterate(data::Mesh, state) = (state += 1; state > lastindex(data) ? nothing : (data[state], state))

Base.values(data::Mesh) = zip(eachcol(data.points), data.measure, data.multiplicity)
Base.firstindex(data::Mesh) = first(CartesianIndices(data.points))[2]
Base.lastindex(data::Mesh) = last(CartesianIndices(data.points))[2]

Base.get(data::Mesh, args...) = Base.getindex(data, args...)
Base.getindex(data::Mesh, i::Int) = (Base.getindex(data.points, :, i), Base.getindex(data.measure, i), Base.getindex(data.multiplicity, i))
Base.getindex(data::Mesh, range::AbstractRange) = (getindex(data, i) for i in range)

points(ks::AbstractMatrix) = ks
points(mesh::Mesh) = mesh.points
points(ks::AbstractVector{<:AbstractVector}) = hcat(ks...)
points(ks::Paths.DiscretePath) = ks.points

################################################################################
################################################################################

function regulargrid(; nk::Int = 100, dim::Int = 2)
    nk = floor(Int, nk^(1 / dim))

    it1d = range(0.0, 1.0; length = nk + 1)[1:end-1]
    itNd = Iterators.product(Iterators.repeated(it1d, dim)...)

    out = hcat(([v...] for v in itNd)...)

    out
end

function randomgrid(; nk::Int = 100, dim::Int = 2, rot_symmetry::Int = 1, B = :id)
    @assert rot_symmetry > 0

    if B == :id
        B = Matrix(1.0 * I, dim, dim)
    end

    N = div(nk, rot_symmetry) # integer division
    @assert N > 0

    println("Random (symmetrized) k grid. # points: ", N * rot_symmetry)
    ks = rand(Float64, (dim, N))

    # symmetrized sampling
    for i = 1:(rot_symmetry-1)
        ks = Matrix([ks inv(B[1:dim, 1:dim]) * rotation2D(2π / rot_symmetry * i) * B[1:dim, 1:dim] * ks])
    end

    ks
end

################################################################################
# Possibly for the future:
################################################################################

# mutable struct MeshData{A<:AbstractMatrix,B<:AbstractVector,C<:AbstractVector{<:Int},D<:AbstractMatrix}
#     points::A # each column is a point
#     measure::B # each point has a measure of it's corresponding plaquette
#     values::D # each point can have multiple values attached (in our case energy bands)
#     multiplicity::C # in case we have to keep track of a hidden multiplicity
# end

# meshweights(data::MeshData) = data.measure .* data.multiplicity

# Base.length(data::MeshData) = Base.size(data.points, 2) # columns in points
# Base.iterate(data::MeshData) = length(data) > 0 ? (data[Base.firstindex(data)], Base.firstindex(data)) : nothing
# Base.iterate(data::MeshData, state) = (state += 1; state > lastindex(data) ? nothing : (data[state], state))

# Base.values(data::MeshData) = zip(eachcol(data.points), eachcol(data.values), data.measure, data.multiplicity)
# Base.firstindex(data::MeshData) = first(CartesianIndices(data.points))[2]
# Base.lastindex(data::MeshData) = last(CartesianIndices(data.points))[2]

# Base.get(data::MeshData, args...) = Base.getindex(data, args...)
# Base.getindex(data::MeshData, i::Int) = (Base.getindex(data.points, :, i), Base.getindex(data.values, :, i), Base.getindex(data.measure, i), Base.getindex(data.multiplicity, i))
# Base.getindex(data::MeshData, range::AbstractRange) = (getindex(data, i) for i in range)

# function getfilling(banddata::MeshData, μ::T1; T::T1 = 0.01) where {T1<:Real}

#     Statistics.mean(sum(m * dk * fermidirac(e - μ; T = T) for (_, e, dk, m) in values(banddata)))
# end