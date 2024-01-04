
import LatticeQM.Operators 
import LatticeQM.TightBinding
# using ..Operators: densitymatrix, densitymatrix!
# using ..TightBinding: zerokey
import SparseArrays

function Operators.getdensitymatrix!(ρ::BdGOperator{T}, H::BdGOperator{T}, ks::AbstractMatrix, μ::Real=0; kwargs...) where {T}
    d = hopdim(H)

    Operators.addchemicalpotential!(H, -μ)
    out = Operators.getdensitymatrix!(ρ.h, H.h, ks, 0.0; kwargs...)
    Operators.addchemicalpotential!(H, μ)
    # TODO: note that out might pick up an additional term - mu * n_el due to the chemical-potential shift
    #       ... figure out if this is a problem or not 

    # MU = Hops(Dict(TightBinding.zerokey(H) => μ .* complex.(SparseArrays.sparse(I,d,d))))
    # addelectronsector!(H, -1 * MU)
    # out = Operators.getdensitymatrix!(ρ.h, H.h, ks, 0.0; kwargs...)
    # addelectronsector!(H, MU)

    # H[TightBinding.zerokey(H)][d+1:2*d, d+1:2*d] .+= 2 * μ * I
    # out = Operators.getdensitymatrix!(ρ.h, H.h, ks, μ; kwargs...)
    # H[TightBinding.zerokey(H)][d+1:2*d, d+1:2*d] .-= 2 * μ * I

    out
end