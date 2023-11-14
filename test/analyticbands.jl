module analyticbands
using LinearAlgebra: dot
export grapheneconductionband

function grapheneconductionband(k)
    rs = [[1, 0], [-1, 0], 
          [0, 1], [0, -1], 
          [1, -1], [-1, 1]]
    total = 3.0
    for r in rs
        term = cos(2 * π * dot(k, r))
        total += term
    end
    energy = sqrt(total)
    return energy
end

end