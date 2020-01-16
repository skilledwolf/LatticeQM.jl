using ..Utils: regulargrid

solvehartreefockmf_k(v::AnyHops, hops::AnyHops, args...; kwargs...) = solvehartreefockmf_k(v, getbloch(hops), args...; kwargs...)
function solvehartreefockmf_k(v::AnyHops, h::Function, ρ_init::AnyHops, filling::Number, args...; kwargs...)
    ℋ_op, ℋ_scalar = hartreefock(h, v)

    solveselfconsistent(ρ_init, ℋ_op, ℋ_scalar, filling, args...; kwargs...)
end

solveselfconsistent(hf, ρ_init::AnyHops, filling::Number, args...; kwargs...) = solveselfconsistent(ρ_init, hf..., ks, filling; kwargs...)

solveselfconsistent(ρ_init, args...; kwargs...) = solveselfconsistent!(copy(ρ_init), args...; kwargs...)
solveselfconsistent!(ρ0::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, args...; kwargs...) = solveselfconsistent!(ρ0, copy(ρ0), ℋ_op, ℋ_scalar, filling, args...; kwargs...)
function solveselfconsistent!(ρ0::AnyHops, ρ1::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64; klin::Int, kwargs...)
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling, regulargrid(nk=klin^2); kwargs...)
end

# abstract type MFobject end
#
# mutable struct HartreeFockMeanField <: MFobject
#     ρ::AnyHops{AbstractMatrix}
#     μ::Float64
# end
#
# SharedArrays.SharedArray(d::AnyHops) = Hops(a=>SharedArray(b) for (a, b)=d)
# SharedArrays.SharedArray(o::HartreeFockMeanField) = SharedArray(o.ρ)

function solveselfconsistent!(ρ0::AnyHops, ρ1::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, ks::AbstractMatrix{Float64};
    iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...)
    """
        Searches a self-consistent meanfield solution for the functional

            ℋ: ρ → h  where h(k) is a hermitian N × N Matrix

        at given filling (between 0 and 1). k space is discretized with the given points ks.

        returns (1) the density matrix of the meanfield (2) ground state energy of the meanfield operator (3) the chemical potential (4) convergence flag (bool) (5) error estimate
    """
    ρ0 = Dict(δL=>SharedArray(m) for (δL, m)=ρ0)
    ρ1 = Dict(δL=>SharedArray(m) for (δL, m)=ρ1)

    function update!(ρ1, ρ0)
        # Update meanfield Hamiltonian and chemical potential
        h = ℋ_op(ρ0)
        Σ = spectrum(h; format=:dense) # lazy diagonalization

        if verbose
            @info("Updating chemical potential for given filling.")
        end
        μ = chemicalpotential(h, ks, filling; T=T)#; format=format) # @time

        # Obtain the meanfield density matrix of the updated Hamiltonian
        if verbose
            @info("Updating the meanfield density matrix.")
        end
        ϵ0 = densitymatrix!(ρ1, Σ, ks, μ; T=T) # @time

        ϵ0 # return the groundstate energy (density matrix was written to ρ1)
    end

    # Compute the ground state energy for the mean-field fixed point
    ϵ0, error, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, kwargs...)

    h = ℋ_op(ρ1)
    μ = chemicalpotential(h, ks, filling; T=T) # Calculate the chemical potential at the end of iteration

    ϵ_offset = ℋ_scalar(ρ1)
    ϵ_GS = ϵ0 + ϵ_offset

    ρ1, ϵ_GS, μ, converged, error
end
# function solveselfconsistent!(ρ0::AnyHops, ρ1::AnyHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, ks::AbstractMatrix{Float64};
#     iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...)
#     """
#         Searches a self-consistent meanfield solution for the functional
#
#             ℋ: ρ → h  where h(k) is a hermitian N × N Matrix
#
#         at given filling (between 0 and 1). k space is discretized with the given points ks.
#
#         returns (1) the density matrix of the meanfield (2) ground state energy of the meanfield operator (3) the chemical potential (4) convergence flag (bool) (5) error estimate
#     """
#     ρ0 = Dict(δL=>SharedArray(m) for (δL, m)=ρ0)
#     ρ1 = Dict(δL=>SharedArray(m) for (δL, m)=ρ1)
#
#     function update!(ρ1, ρ0)
#         # Update meanfield Hamiltonian and chemical potential
#         h = ℋ_op(ρ0)
#         Σ = spectrum(h; format=:dense) # lazy diagonalization
#
#         if verbose
#             @info("Updating chemical potential for given filling.")
#         end
#         μ = chemicalpotential(h, ks, filling; T=T)#; format=format) # @time
#
#         # Obtain the meanfield density matrix of the updated Hamiltonian
#         if verbose
#             @info("Updating the meanfield density matrix.")
#         end
#         ϵ0 = densitymatrix!(ρ1, Σ, ks, μ; T=T) # @time
#
#         ϵ0 # return the groundstate energy (density matrix was written to ρ1)
#     end
#
#     # Compute the ground state energy for the mean-field fixed point
#     ϵ0, error, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, kwargs...)
#
#     h = ℋ_op(ρ1)
#     μ = chemicalpotential(h, ks, filling; T=T) # Calculate the chemical potential at the end of iteration
#
#     ϵ_offset = ℋ_scalar(ρ1)
#     ϵ_GS = ϵ0 + ϵ_offset
#
#     ρ1, ϵ_GS, μ, converged, error
# end

###################################################################################################
# Backwards compatibility
###################################################################################################
export solve_selfconsistent
@legacyalias solveselfconsistent solve_selfconsistent
