import LatticeQM.Meanfield
# using ..Meanfield: hartreefock_pairing
import LatticeQM.Utils

mutable struct HartreeFockBDG{K,T2,T<:Hops{K,T2}} <: Meanfield.MeanfieldGenerator{T}
    h::T
    v::T
    μ::Float64
    V0::T2 # Assuming T2 is Float64, adjust accordingly
    hMF::T
    ΔMF::T
    ϵMF::Float64
    fock::Bool
    hartree::Bool

    function HartreeFockBDG(h::T, v::T, μ=0.0; hartree=true, fock=true) where {K,T2,T<:Hops{K,T2}}
        # Extract K and T2 from T if needed inside the constructor
        # Example: Suppose AbstractHops is defined as AbstractHops{K, T2}
        # K = Base.parameter_upper_bound(T, 1)
        # T2 = Base.parameter_upper_bound(T, 2)
        V0 = sum(v[L] for L in keys(v))
        hMF = zero(h)  # Pre-allocated for performance
        ΔMF = zero(h)  # Pre-allocated for performance
        ϵMF = 0.0
        new{K,T2,T}(h, v, μ, V0, hMF, ΔMF, ϵMF, fock, hartree)
    end
end

function Meanfield.HartreeFock(h::BdGOperator{T}, v::T, μ=0.0; hartree=true, fock=true) where {K,T2,T<:Hops{K,T2}}
    h_el = Utils.getelectronsector(h)
    HartreeFockBDG(h_el, v, μ; hartree=hartree, fock=fock)
end

Meanfield.hMF(hf::HartreeFockBDG) = BdGOperator(hf.hMF, hf.ΔMF)

function (hf::HartreeFockBDG)(ρ) # where {K,T2,T<:Hops{K,T2}} #(ρ::T) where {K,T2,T<:Hops{K,T2}}
    ρ_el = getelectronview(ρ)
    ρ_Δ = getpairingview(ρ)
    meanfieldOperator!(hf, ρ_el, ρ_Δ)
    meanfieldScalar!(hf, ρ_el, ρ_Δ)
    hf
end

function initialize_ΔMF!(hf::HartreeFockBDG, ρΔ)
    # Initialize hf.ΔMF
    # Adjust this part according to the actual structure and requirements of hf.ΔMF
    T = eltype(hf.h[zerokey(hf.h)])
    dim = size(hf.h[zerokey(hf.h)])
    for k in union(keys(hf.h), keys(hf.v), keys(ρΔ)) # Adjust this line based on how hf.ΔMF should be initialized
        if haskey(hf.ΔMF, k)
            hf.ΔMF[k] .= 0.0 # Assuming default initialization to 0.0
        else
            hf.ΔMF[k] = zeros(T, dim) # Adjust this line as needed
        end
    end
end

function meanfieldOperator!(hf::HartreeFockBDG, ρ_el, ρΔ)
    Meanfield.initialize_hMF!(hf, ρ_el)
    initialize_ΔMF!(hf, ρΔ)

    if hf.fock
        Meanfield.meanfieldOperator_addfock!(hf, ρ_el)
        Meanfield.meanfieldOperator_addfock_pairing!(hf, ρΔ)
    end
    if hf.hartree
        Meanfield.meanfieldOperator_addhartree!(hf, ρ_el)
    end
    nothing
end

function meanfieldScalar!(hf::HartreeFockBDG, ρ_el, ρΔ)
    hf.ϵMF = 0.0

    if hf.fock
        hf.ϵMF += Meanfield.meanfieldScalar_fock(hf, ρ_el)
        hf.ϵMF += Meanfield.meanfieldScalar_fock_pairing(hf, ρΔ)
    end
    if hf.hartree
        hf.ϵMF += Meanfield.meanfieldScalar_hartree(hf, ρ_el)
    end
    nothing
end

