### Utility functions

α(n::Int; m::Int=1) = acos((3.0*n^2 + 3*n*m + m^2/2.0)/(3.0*n^2 + 3*n*m + m^2))

# using .Structure: translate!, fold_atoms!, mirrorZ!, build_superlattice, add_dimension!, repeat_atoms!, rotate_atoms_XY!, crop_to_unitcell!, combine_lattices!

function twist_triangular_2D(triang_lat::Lattice{T}, triang_lat2::Lattice{T}, n::Int; z::Float64=3, m::Int=1) where {T<:AbstractMatrix{Float64}}
    """
    This function assumes that the triangular 2D lattice "triang_lat1,2" has at least one layer
    with z=0 (and possibly more layers with z>0).
    triang_lat1 and triang_lat2 must have identical layers at z=0
    """
    triang_lat2 = deepcopy(triang_lat2) # just to be save

    # Crucual parameters
    twist_angle = α(n)
    superperiods = [[n; n+m] [-n-m; 2*n+m]]

    # Move the initial layer up along z away from z=0
    translate!(triang_lat, "z", z/2)

    # Copy the layer and start from AB stacking at the twist interface
    # (the construction demands it for some reason)
    # Then mirror the layer at z=0 plane
    translate!(triang_lat2, "z", z/2)
    triang_lat2.atoms .= - triang_lat2.atoms
    fold_atoms!(triang_lat2)
    mirrorZ!(triang_lat2)

    # Build (non-orthogonal) supercells and move it up along z
    superlat1 = build_superlattice(triang_lat, superperiods)
    add_dimension!(superlat1, "layer", fill(0.0, (1, atom_count(superlat1))))
    superlat2 = build_superlattice(triang_lat2, superperiods)
    add_dimension!(superlat2, "layer", fill(1.0, (1, atom_count(superlat2))))

    # Make sure all points lie in the primitive unit cell
    fold_atoms!(superlat1)
    fold_atoms!(superlat2)

    # Rotate the atom positions (keeping the lattice vectors fixed)
    @info "Twist angle α="*string(round(twist_angle/π*180; digits=3))*"°"
    repeat_atoms!(superlat2, [-1:1,-1:1])
    rotate_atoms_XY!(superlat2, twist_angle)
    crop_to_unitcell!(superlat2)

    combine_lattices!(superlat1, superlat2)

    superlat1
end
