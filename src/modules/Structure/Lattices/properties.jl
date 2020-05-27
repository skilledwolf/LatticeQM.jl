###################################################################################################
# Properties
###################################################################################################

latticedim(lat::Lattice) = lat.latticedim
countorbitals(lat::Lattice) = size(lat.spacecoordinates, 2)
extraspacedim(lat::Lattice) = size(lat.extracoordinates, 1)
spacedim(lat::Lattice) = size(getA(lat), 1) # used to be latticedim(lat) + extraspacedim(lat)

hasdimension(lat::Lattice, name::String) = haskey(lat.extralabels, name)
assertdimension(lat::Lattice, name::String) = !hasdimension(lat, name) ? error("No $name coordinates specified.") : Nothing

basis(lat::Lattice, selector=(:)) = lat.basis[:,selector]
getA(lat::Lattice, selector=(:)) = basis(lat)[:,1:latticedim(lat)][:,selector]


"""
Calculate the dual lattice of lat.A.
Here we use the general formula \$B = A * (A^T * A)^-1\$. That works also 
when the d-dim lattice is embedded in D-dim space.
"""
function getB(lat::Lattice, selector=(:))
    d = latticedim(lat)
    A = getA(lat)[:,1:d]
    return (A * inv(transpose(A)*A))[:,selector]
end

coordinates(lat::Lattice, selector=(:)) = lat.spacecoordinates[selector,:]
# positions(lat::Lattice) = getA(lat) * coordinates(lat)
positions(lat::Lattice, selector=(:)) = (basis(lat) * coordinates(lat))[selector,:]


allpositions(lat::Lattice, args...) = [ positions(lat); extracoordinates(lat, args...)]

extracoordinates(lat::Lattice) = lat.extracoordinates
extracoordinates(lat::Lattice, dim::String) = extracoordinates(lat, [dim])
function extracoordinates(lat::Lattice, dims::Vector{String})
    for dim in dims
        assertdimension(lat, dim)
    end
    indices = [lat.extralabels[dim] for dim=dims]
    extracoordinates(lat)[indices, :]
end

function setextracoordinates!(lat::Lattice, dim::String, values)
    assertdimension(lat, dim)
    lat.extracoordinates[lat.extralabels[dim],:] = values
end

function filterindices(lat::Lattice, name::String, condition::Function)
    return [index for (index, val) in enumerate(extracoordinates(lat, name)) if condition(val)]
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


function Base.show(io::IO, lat::Lattice)
    println(io, "Lattice dimension:     ", latticedim(lat))
    println(io, "Space dimension:       ", spacedim(lat))
    println(io, "Number of atoms:       ", countorbitals(lat))
    println(io, "Non-spatial dimension: ", extraspacedim(lat))
    println(io, "Basis:")
    show(io, basis(lat))
    println(io, "")
    println(io, "Orbital/atom coordinates: ")
    show(io, allpositions(lat))
end

function Base.show(io::IO, m::MIME"text/plain", lat::Lattice)
    println(io, "Lattice dimension:     ", latticedim(lat))
    println(io, "Space dimension:       ", spacedim(lat))
    println(io, "Number of atoms:       ", countorbitals(lat))
    println(io, "Non-spatial dimension: ", extraspacedim(lat))
    println(io, "Basis:")
    show(io, m, basis(lat))
    println(io, "")
    println(io, "Orbital/atom coordinates: ")
    show(io, m, allpositions(lat))
end
