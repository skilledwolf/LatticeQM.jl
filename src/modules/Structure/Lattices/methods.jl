@legacyalias foldcoordinates! fold_atoms!
function foldcoordinates!(lat::Lattice)
    lat.orbitalcoordinates[:] .= (mod.(lat.orbitalcoordinates, 1))[:]
end

@legacyalias rotatebasis! rotate_A_XY!
function rotatebasis!(lat::Lattice, α::Float64)
    @assert latticedim(lat) == 2
    lat.latticevectors[1:2,1:2] .= (RotationMatrix(α) * getA(lat)[1:2,1:2])
end

@legacyalias rotatecoordinates! rotate_atoms_XY!
function rotatecoordinates!(lat::Lattice, θ::Float64)
    mat = inv(getA(lat)) * RotationMatrix(θ)* getA(lat)
    lat.orbitalcoordinates[1:2,:] .= mat * lat.orbitalcoordinates[1:2,:]
end

function translate!(lat::Lattice, name::String, δ::Float64)
    @assert hasdimension(lat, name)
    lat.extrapositions[lat.extradimensions[name],:] .+= δ
    lat
end

mirrorZ(lat::Lattice) = mirrorZ!(copy(lat))
function mirrorZ!(lat::Lattice)
    @assert hasdimension(lat, "z")
    lat.extrapositions[lat.extradimensions["z"],:] .*= (-1.0)
    lat
end

@legacyalias newdimension! add_dimension!
newdimension!(lat::Lattice, name::String, extrapositions::AbstractVector{Float64}) = newdimension!(lat,name,extrapositions')
function newdimension!(lat::Lattice, name::String, extrapositions::AbstractMatrix{Float64})
    @assert size(lat.orbitalcoordinates)[2] == size(extrapositions)[2]

    lat.extradimensions[name] = extraspacedim(lat)+1
    lat.extrapositions = vcat(lat.extrapositions, extrapositions)
    lat
end

@legacyalias mergelattices combine_lattices
@legacyalias mergelattices! combine_lattices!
mergelattices(lat1::Lattice, lat2::Lattice) = mergelattices!(copy(lat1),lat2)
function mergelattices!(lat1::Lattice, lat2::Lattice)
    @assert getA(lat1) ≈ getA(lat2)
    @assert lat1.extradimensions == lat2.extradimensions
    lat1.orbitalcoordinates = hcat(lat1.orbitalcoordinates, lat2.orbitalcoordinates)
    lat1.extrapositions = hcat(lat1.extrapositions, lat2.extrapositions)

    lat1
end


