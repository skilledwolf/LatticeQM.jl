###################################################################################################
# Properties
###################################################################################################

latticedim(lat::Lattice) = size(getA(lat), 2)
countorbitals(lat::Lattice) = size(lat.orbitalcoordinates, 2)
extraspacedim(lat::Lattice) = size(lat.extrapositions, 1) - latticedim(lat) # used to be size(lat.extrapositions, 1)
spacedim(lat::Lattice) = size(getA(lat), 1) # used to be latticedim(lat) + extraspacedim(lat)

hasdimension(lat::Lattice, name::String) = haskey(lat.extradimensions, name)
assertdimension(lat::Lattice, name::String) = !hasdimension(lat, name) ? error("No $name coordinates specified.") : Nothing

getA(lat::Lattice) = lat.latticevectors

"""
Calculate the dual lattice of lat.A.
Here we use the general formula \$B = A * (A^T * A)^-1\$. That works also 
when the d-dim lattice is embedded in D-dim space.
"""
function getB(lat::Lattice)
    A = getA(lat)
    return A * inv(transpose(A)*A)
end

coordinates(lat::Lattice) = lat.orbitalcoordinates
# positions(lat::Lattice) = getA(lat) * coordinates(lat)
function positions(lat::Lattice)
    d = latticedim(lat)
    D = spacedim(lat)
    A = getA(lat)

    X = coordinates(lat)
    R = extrapositions(lat)[1:D-d,:]

    R0 = A * X
    R0[d+1:D,:] .+= R

    R0
end

function allpositions(lat::Lattice, args...)
    d = latticedim(lat)
    D = spacedim(lat)

    [ positions(lat); extrapositions(lat, args...)[1+D-d:end,:] ]
end

extrapositions(lat::Lattice) = lat.extrapositions#[spacedim(lat)+1:end,:]
extrapositions(lat::Lattice, dim::String) = extrapositions(lat, [dim])
function extrapositions(lat::Lattice, dims::Vector{String})
    for dim in dims
        assertdimension(lat, dim)
    end
    indices = [lat.extradimensions[dim] for dim=dims]
    extrapositions(lat)[indices, :]
end

function setextrapositions!(lat::Lattice, dim::String, values)
    assertdimension(lat, dim)
    lat.extrapositions[lat.extradimensions[dim],:] = values
end

function filterindices(lat::Lattice, name::String, condition::Function)
    return [index for (index, val) in enumerate(extrapositions(lat, name)) if condition(val)]
end
filterpositions(lat::Lattice, args...) = positions(lat)[:,filterindices(lat, args...)]
filtercoordinates(lat::Lattice, args...) = coordinates(lat)[:,filterindices(lat, args...)]

function fractionalize(lat::Lattice, positions::AbstractMatrix{Float64})
    A = getA(lat)
    return transpose(A * inv(transpose(A) * A)) * positions
end
function fractionalize!(lat::Lattice, positions::AbstractMatrix{Float64})
    A = getA(lat)
    d = latticedim(lat)

    positions[1:d,:] .= (transpose(A * inv(transpose(A) * A)) * positions)
end
foldfractional(fracpositions::AbstractMatrix{Float64}) = mod.(fracpositions,1) 

