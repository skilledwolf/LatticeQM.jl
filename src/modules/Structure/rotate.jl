import LinearAlgebra: norm, dot, cross

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

    # Rodrigues form, https://en.wikipedia.org/wiki/Rotation_matrix
    # NOTE: a rotation matrix is NOT symmetric — the off-diagonal partners
    # differ by ±2uᵢs. The old code mirrored r12/r13/r23 across the diagonal,
    # producing a non-orthogonal matrix (det = cos2θ for n = ẑ) that shrank
    # generic vectors and broke the spiral initial guess.
    c = cos(θ); s = sin(θ)
    d1 = c + u1^2 * (1-c)
    d2 = c + u2^2 * (1-c)
    d3 = c + u3^2 * (1-c)

    [ d1                    u1*u2*(1-c) - u3*s    u1*u3*(1-c) + u2*s;
      u2*u1*(1-c) + u3*s    d2                    u2*u3*(1-c) - u1*s;
      u3*u1*(1-c) - u2*s    u3*u2*(1-c) + u1*s    d3 ]
end


"""
    signedangle(e1::T,e2::T; z=nothing)

returns the signed angle between 3D vectors e1 and e2.
z should be the normal vector.

"""
function signedangle(e1::T,e2::T; z=[0.0,0,1]) where T<:AbstractVector
    @assert length(e1)==length(e2)==3 "Vectors must have length 3."
    e1 = e1/norm(e1)
    e2 = e2/norm(e2)
    z = z/norm(z)

    atan(dot(cross(e1,e2),z), dot(e1,e2))
end