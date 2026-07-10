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

σUP = (σ0+σ3)/2
σDOWN = (σ0-σ3)/2
σPLUS = (σ1+1im*σ2)/2
σMINUS = (σ1-1im*σ2)/2

# σs = [σ1, σ2, σ3]
σs = Dict(
      0=>σ0, 1=> σ1, 2=> σ2, 3=> σ3,
      "X"=>σ1, "Y"=>σ2, "Z"=>σ3,
      "up"=>σUP, "down"=>σDOWN,
      "plus"=>σPLUS, "minus"=>σMINUS
)

import LinearAlgebra: norm

"""
    spinorrotation(θ, n=[0,0,1])

Returns SU(2) spinor rotation matrix U parametrized by rotation angle θ and rotation axis n=[n1,n2,n3].
U = exp(-i (θ/2) σ.n).
    
Note: This routine will always normalize n.
"""
function spinorrotation(θ, n::Vector=[0,0,1])
      @assert length(n) == 3 "Normal vector n for rotation must have 3 components."
      n̂ = n / norm(n)
      M = n̂[1]*σ1 + n̂[2]*σ2 + n̂[3]*σ3
      σ0 * cos(θ/2) - 1im * M * sin(θ/2)
end