using ..Structure: regulargrid

# Specialized cases
solvehartreefock(h, v, ρ_init, filling::Number, args...; kwargs...) = solveselfconsistent(ρ_init, hartreefock(h, v)..., filling, args...; kwargs...)
solvefock(h, v, ρ_init, filling::Number, args...; kwargs...) = solveselfconsistent(ρ_init, fock(h, v)..., filling, args...; kwargs...)
solvehartree(h, v, ρ_init, filling::Number, args...; kwargs...) = solveselfconsistent(ρ_init, hartree(h, v), filling, args...; kwargs...)

# Interface to solveselfconsistent!(ρ0, ...)
solveselfconsistent(ρ0, ℋ_op, ℋ_scalar, filling::Number, ks; kwargs...) = solveselfconsistent!(deepcopy(ρ0), ℋ_op, ℋ_scalar, filling, ks; kwargs...)
solveselfconsistent(ρ0, ℋ_op, ℋ_scalar, filling::Number; klin, kwargs...) = solveselfconsistent!(deepcopy(ρ0), ℋ_op, ℋ_scalar, filling; klin=klin, kwargs...)

# Interface to solveselfconsistent!(ρ0, ρ1, ...)
function solveselfconsistent!(ρ0, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, args...; kwargs...)
    return solveselfconsistent!(ρ0, deepcopy(ρ0), ℋ_op, ℋ_scalar, filling, args...; kwargs...)
end

# Interface to solveselfconsistent!(ρ0, ρ1, ..., ks, ...)
function solveselfconsistent!(ρ0, ρ1, ℋ_op::Function, ℋ_scalar::Function, filling::Float64; klin::Int, kwargs...)
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling, regulargrid(nk=klin^2); kwargs...)
end

mutable struct Hamiltonian{T}
    h::T
    μ::Float64
end

# using JLD
import Distributed: nprocs
import ..TightBinding: SharedDenseHops
import ..Operators: getdensitymatrix!

"""
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling, ks; convergenceerror=false, multimode=:serial, checkpoint::String="", hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...)
    solveselfconsistent!(ρ0, ℋ_op, ℋ_scalar, filling, ks; kwargs...)
    solveselfconsistent!(ρ0, ρ1, ℋ_op, ℋ_scalar, filling; klin, kwargs...)
    solveselfconsistent!(ρ0, ℋ_op, ℋ_scalar, filling; klin, kwargs...)

Searches a self-consistent meanfield solution for the functional ℋ: ρ → h
at given filling (between 0 and 1). k space is discretized with the given points ks.
returns (1) the density matrix of the meanfield (2) ground state energy of the meanfield operator (3) the chemical potential (4) convergence flag (bool) (5) error estimate

If the checkpoint keyword is set, e.g. `checkpoint="mf.jld"`, the mean field will be saved into a file at each step of the iteration,
allowing to interrupt and restart a long-running calculation safely. If `checkpoint` is set and the file exists, this method
will assume that the file contains a valid mean-field and use it as initial guess. To avoid this behavior, specify `hotstart=false`.

parallel=true might help if diagonalization per k point is very time consuming
(e.g. for twisted bilayer graphene)
note that for small problems `parallel=true` may decrease performance (communication overhead)
"""
function solveselfconsistent!(ρ0::AbstractHops, ρ1::AbstractHops, ℋ_op::Function, ℋ_scalar::Function, filling::Float64, ks;
    convergenceerror=false, multimode=:serial, checkpoint::String="", callback=(x->nothing), hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...)

    # if checkpoint != "" && isfile(checkpoint) && hotstart
    #     println("Loading checkpoint file as initial guess: $checkpoint")
    #     ρ0 = JLD.load(checkpoint, "mf")
    # end

    # Turn dense and prepare for distributed computing
    ρ0 = (multimode==:distributed && nprocs()>1) ? SharedDenseHops(ρ0) : DenseHops(ρ0)
    ρ1 = (multimode==:distributed && nprocs()>1) ? SharedDenseHops(ρ1) : DenseHops(ρ1)

    H = Hamiltonian(ℋ_op(ρ0), 0.0)

    function updateH!(H::Hamiltonian, ρ)
        verbose ? @info("Updating chemical potential for given filling...") : nothing
        
        H.h = ℋ_op(ρ) # get updated Hamiltonian
        H.μ = chemicalpotential(H.h, ks, filling; T=T, multimode=multimode)
        H
    end

    function update!(ρ1, ρ0)
        # @time updateH!(H, ρ0)
        updateH!(H, ρ0)

        verbose ? @info("Updating the mean field...") : nothing
        # @time ϵ0 = getdensitymatrix!(ρ1, H.h, ks, H.μ; multimode=multimode, T=T, format=:dense) # get new meanfield and return the groundstate energy (density matrix was written to ρ1)
        ϵ0 = getdensitymatrix!(ρ1, H.h, ks, H.μ; multimode=multimode, T=T, format=:dense) # get new meanfield and return the groundstate energy (density matrix was written to ρ1)

        callback(ρ1)

        # if checkpoint != ""
        #     verbose ? @info("Saving intermediate mean field...") : nothing
        #     JLD.save(checkpoint, "mf", DenseHops(ρ1))
        # end

        ϵ0 + ℋ_scalar(ρ1) # proper ground state energy
    end

    # Compute the ground state energy for the mean-field fixed point
    ϵ_GS, residual, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, verbose=verbose, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end
    updateH!(H, ρ1)

    DenseHops(ρ1), ϵ_GS, H, converged, residual
end

