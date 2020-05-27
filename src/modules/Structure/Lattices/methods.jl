function foldcoordinates!(lat::Lattice)
    d = latticedim(lat)
    lat.spacecoordinates[1:d,:] = mod.(lat.spacecoordinates[1:d,:], 1)
end

function rotatebasis!(lat::Lattice, α::Float64)
    @assert spacedim(lat) >= 2
    lat.basis[1:2,:] .= (rotation2D(α) * basis(lat)[1:2,:])
end

function rotatecoordinates!(lat::Lattice, θ::Float64)
    R = positions(lat); D = spacedim(lat); d = latticedim(lat)
    R[1:2,:] .= rotation2D(θ) * R[1:2,:]
    fractionalize!(lat, R)

    lat.spacecoordinates[:,:] .= R[:,:]
    # lat.extracoordinates[1:D-d,:] .= R[d+1:D,:]
end

function translate!(lat::Lattice, name::String, δ::Float64)
    @assert hasdimension(lat, name)
    lat.extracoordinates[lat.extralabels[name],:] .+= δ
    lat
end

function translate!(lat::Lattice, i::Integer, δ::Float64)
    lat.spacecoordinates[:,:] .+= basis(lat,i) * δ
    lat
end

function translate!(lat::Lattice, v::AbstractVector)
    lat.spacecoordinates[:,:] .+= basis(lat) * v
    lat
end

"""
    displace!(lat, f::Function)

Function f takes orbital i at position p_i and displaces it by vector v_i = f(p_i).
"""
function displace!(lat, f::Function)

    # Modify the lattice geometry
    XY = positions(lat)

    for (i_,p)=enumerate(eachcol(XY))
        lat.spacecoordinates[:,i_] += f(p)
    end

    lat
end

"""
    displaceZ!(lat, f::Function)

Function f takes orbital i at position p_i and displaces it in the third coodinate by z_i = f(p_i).
"""
function displaceZ!(lat, f::Function)

    # Modify the lattice geometry
    XY = positions(lat)

    for (i_,p)=enumerate(eachcol(XY))
        lat.spacecoordinates[3,i_] += f(p)
    end

    lat
end

mirrorZ(lat::Lattice) = mirrorZ!(copy(lat))
function mirrorZ!(lat::Lattice)
    lat.spacecoordinates[3,:] .*= (-1.0)
    lat
end
# function mirrorZ!(lat::Lattice)
#     @assert hasdimension(lat, "z")
#     lat.extracoordinates[lat.extralabels["z"],:] .*= (-1.0)
#     lat
# end

newdimension!(lat::Lattice, name::String, extracoordinates::AbstractVector{Float64}) = newdimension!(lat,name,extracoordinates')
function newdimension!(lat::Lattice, name::String, extracoordinates::AbstractMatrix{Float64})
    @assert countorbitals(lat) == size(extracoordinates,2)

    lat.extralabels[name] = extraspacedim(lat)+1
    lat.extracoordinates = vcat(lat.extracoordinates, extracoordinates)
    lat
end

mergelattices(lat1::Lattice, lat2::Lattice) = mergelattices!(copy(lat1),lat2)
function mergelattices!(lat1::Lattice, lat2::Lattice)
    @assert getA(lat1) ≈ getA(lat2)
    @assert lat1.extralabels == lat2.extralabels
    lat1.spacecoordinates = hcat(lat1.spacecoordinates, lat2.spacecoordinates)
    lat1.extracoordinates = hcat(lat1.extracoordinates, lat2.extracoordinates)

    lat1
end


