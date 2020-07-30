"""
    rotation2D(θ)

Returns the 2x2 rotation matrix R = R^T parametrized by the rotation angle θ.
"""
rotation2D(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


"""
    rotation3D(θ, n=[0,0,1])

Returns SO(3) rotation matrix parametrized by rotation angle θ and rotation axis n=[n1,n2,n3].
    
Note: This routine will always normalize n.
"""
function rotation3D(θ, n::Vector=[0,0,1])

    u1, u2, u3 = n / norm(n) # normalize n

    # just taken from wikipedia https://en.wikipedia.org/wiki/Rotation_matrix
    c = cos(θ); s = sin(θ)
    d1 = c + u1^2 * (1-c)
    d2 = c + u2^2 * (1-c)
    d3 = c + u3^2 * (1-c)

    r12 = u1 * u2 * (1-c) - u3 * s
    r13 = u1 * u3 * (1-c) + u2 * s
    r23 = u2 * u3 * (1-c) + u1 * s

    [ d1  r12 r13; 
      r12 d2  r23;
      r13 r23 d3 ]
end