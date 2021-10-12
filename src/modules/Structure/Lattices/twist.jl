import LinearAlgebra: dot, norm

function twistangle(n::Int; m::Int=1, degrees=false)
    α = acos((3.0*n^2 + 3*n*m + m^2/2.0)/(3.0*n^2 + 3*n*m + m^2))

    if degrees 
        α = α*180/pi
    end
    
    α
end

twist(lat::Lattice, n::Int; kwargs...) = twist(lat,lat,n; kwargs...)
function twist(lat1::Lattice, lat2::Lattice, n::Int; z::Float64=3, m::Int=1, verbose::Bool=true)
    @assert latticedim(lat1) == 2 && latticedim(lat2) == 2 "twist(...) is only defined for 2D lattices."
    @assert getA(lat1) ≈ getA(lat2) "The two lattices must have the same lattice vectors."
    @assert abs(dot(getA(lat1)[:,1], getA(lat1)[:,2])/(norm(getA(lat1)[:,1])*norm(getA(lat1)[:,1]))) ≈ 0.5 "twist(...) is only defined for triangular lattices."
    """
    This function assumes that the triangular 2D lattice "lat11,2" has at least one layer
    with z=0 (and possibly more layers with z>0).
    lat11 and lat2 must have identical layers at z=0
    """
    lat2 = deepcopy(lat2) # just to be save

    # Crucial parameters
    angle = twistangle(n)
    superperiods = [[n; n+m] [-n-m; 2*n+m]]

    # Move the initial layer up along z away from z=0
    translate!(lat1, 3, z/2) #translate!(lat1, "z", z/2)

    # Copy the layer and start from AB stacking at the twist interface
    # (the construction demands it for some reason)
    # Then mirror the layer at z=0 plane
    translate!(lat2, 3, z/2) # translate!(lat2, "z", z/2)
    lat2.spacecoordinates .*= -1.0
    foldcoordinates!(lat2)
    # mirrorZ!(lat2)

    # Build (non-orthogonal) supercells
    superlat1 = superlattice(lat1, superperiods)
    newdimension!(superlat1, "layer", fill(0.0, (1, countorbitals(superlat1))))
    superlat2 = superlattice(lat2, superperiods)
    newdimension!(superlat2, "layer", fill(1.0, (1, countorbitals(superlat2))))

    # Rotate the atom positions (keeping the lattice vectors fixed)
    if verbose
        println("Twist α="*string(round(angle/π*180; digits=3))*"°   (n,m)=($n,$m)")
    end
    # repeat!(superlat1, [-1:1,-1:1])
    repeat!(superlat2, [-1:1,-1:1])
    # rotatecoordinates!(superlat1, -angle/2)
    rotatecoordinates!(superlat2, angle)
    # crop2unitcell!(superlat1)
    crop2unitcell!(superlat2)

    mergelattices!(superlat1, superlat2)
    foldcoordinates!(superlat1)

    rotatebasis!(superlat1, -angle/2)

    superlat1
end
precompile(twist, (Lattice, Int))
precompile(twist, (Lattice, Lattice, Int))
