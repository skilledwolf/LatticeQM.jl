σ0 = [1.0 0.0;
      0.0 1.0]
σ1 = [0.0 1.0;
      1.0 0.0]
σ2 = [0.0 -1.0im;
      1.0im 0.0]
σ3 = [1.0 0.0;
      0.0 -1.0]

σX = σ1
σY = σ2
σZ = σ3

σs = [σ1, σ2, σ3]

"""
    spinorrotation(θ, n=[0,0,1])

Returns SU(2) spinor rotation matrix U parametrized by rotation angle θ and rotation axis n=[n1,n2,n3].
U = exp(-i (θ/2) σ.n).
    
Note: This routine will always normalize n.
"""
function spinorrotation(θ, n::Vector=[0,0,1])
      @assert length(n) == 3 "Normal vector n for rotation must have 3 components."
      u1, u2, u3 = n / norm(n)

      M = zero(σ2)
      for (u,σ) in zip(n,σs)
            M += (u * σ)
      end

      σ0 * cos(θ/2) - 1im * M * sin(θ/2)
end