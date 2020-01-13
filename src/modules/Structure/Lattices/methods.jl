@legacyalias foldcoordinates! fold_atoms!
function foldcoordinates!(lat::Lattice)
    lat.atoms[:] .= (mod.(lat.atoms, 1))[:]
end

@legacyalias rotatebasis! rotate_A_XY!
function rotatebasis!(lat::Lattice, α::Float64)
    @assert latticedim(lat) == 2
    lat.A[1:2,1:2] .= (RotationMatrix(α) * lat.A[1:2,1:2])
end

@legacyalias rotatecoordinates! rotate_atoms_XY!
function rotatecoordinates!(lat::Lattice, θ::Float64)
    mat = inv(lat.A) * RotationMatrix(θ)* lat.A
    lat.atoms[1:2,:] .= mat * lat.atoms[1:2,:]
end

function translate!(lat::Lattice, name::String, δ::Float64)
    @assert hasdimension(lat, name)
    lat.atoms_aux[lat.extradimensions[name],:] .+= δ
    lat
end

mirrorZ(lat::Lattice) = mirrorZ!(copy(lat))
function mirrorZ!(lat::Lattice)
    @assert hasdimension(lat, "z")
    lat.atoms_aux[lat.extradimensions["z"],:] .*= (-1.0)
    lat
end

@legacyalias newdimension! add_dimension!
newdimension!(lat::Lattice, name::String, extrapositions::AbstractVector{Float64}) = newdimension!(lat,name,extrapositions')
function newdimension!(lat::Lattice, name::String, extrapositions::AbstractMatrix{Float64})
    @assert size(lat.atoms)[2] == size(extrapositions)[2]

    lat.extradimensions[name] = extraspacedim(lat)+1
    lat.atoms_aux = vcat(lat.atoms_aux, extrapositions)
    lat
end

@legacyalias mergelattices combine_lattices
@legacyalias mergelattices! combine_lattices!
mergelattices(lat1::Lattice, lat2::Lattice) = mergelattices!(copy(lat1),lat2)
function mergelattices!(lat1::Lattice, lat2::Lattice)
    @assert lat1.A == lat2.A
    @assert lat1.extradimensions == lat2.extradimensions
    lat1.atoms = hcat(lat1.atoms, lat2.atoms)
    lat1.atoms_aux = hcat(lat1.atoms_aux, lat2.atoms_aux)

    lat1
end


