# using JLD
import LatticeQM.Structure
import LatticeQM.Operators
import LatticeQM.Spectrum
import LatticeQM.TightBinding
using LatticeQM.Utils.Context

# Specialized cases
solvehartreefock(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v), filling, args...; kwargs...)
solvefock(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)

# Interface to solveselfconsistent!(ρ0, ...)
solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling::Number, ks; kwargs...) = solveselfconsistent!(deepcopy(ρ0), mf, filling, ks; kwargs...)
solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling::Number; klin, kwargs...) = solveselfconsistent!(deepcopy(ρ0), mf, filling; klin=klin, kwargs...)

# Interface to solveselfconsistent!(ρ0, ρ1, ...)
solveselfconsistent!(ρ0, mf::MeanfieldGenerator, filling::Float64, args...; kwargs...) = solveselfconsistent!(ρ0, deepcopy(ρ0), mf, filling, args...; kwargs...)

# Interface to translate klin to solveselfconsistent!(ρ0, ρ1, ..., ks, ...)
solveselfconsistent!(ρ0, ρ1, mf::MeanfieldGenerator, filling::Float64; klin::Int, kwargs...) = solveselfconsistent!(ρ0, ρ1, mf, filling, Structure.regulargrid(nk=klin^2); kwargs...)

function solveselfconsistent!(ρ0::T, ρ1::T, mf::MeanfieldGenerator, filling::Float64, ks; multimode=:global, kwargs...) where {T}
    c = trycontext(ensurecontext(multimode), SerialContext())
    solveselfconsistent!(c, ρ0, ρ1, mf::MeanfieldGenerator, filling::Float64, ks; kwargs...)
end

solveselfconsistent!(::SerialContext, ρ0::T, ρ1::T, args...; kwargs...) where {T} = solveselfconsistent!(DummyContext(), TightBinding.dense(ρ0), TightBinding.dense(ρ1), args...; multimode=:serial, kwargs...)
solveselfconsistent!(::DistributedContext, ρ0::T, ρ1::T, args...; kwargs...) where {T} = solveselfconsistent!(DummyContext(), TightBinding.shareddense(ρ0), TightBinding.shareddense(ρ1), args...; multimode=:distributed, kwargs...)

function solveselfconsistent!(::MultiThreadedContext, args...; kwargs...)
    error("Multithreaded context not implemented yet.")
end

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
function solveselfconsistent!(::DummyContext, ρ0::T1, ρ1::T1, hartreefock::MeanfieldGenerator, filling::Float64, ks;
    convergenceerror=false, multimode=:serial, checkpoint::String="", callback=(x -> nothing), hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false, kwargs...) where {T1}

    # if checkpoint != "" && isfile(checkpoint) && hotstart
    #     println("Loading checkpoint file as initial guess: $checkpoint")
    #     ρ0 = JLD.load(checkpoint, "mf")
    # end

    # println("rho0: ", typeof(ρ0))
    # println("rho1: ", typeof(ρ1))

    function update!(ρ1, ρ0)
        
        hartreefock(ρ0) # update meanfield (h is updated in-place)

        verbose ? @info("Updating chemical potential for given filling...") : nothing
        hartreefock.μ = Spectrum.chemicalpotential(hMF(hartreefock), ks, filling; T=T, multimode=multimode)


        verbose ? @info("Updating the mean field density matrix...") : nothing
        ϵ0 = Operators.getdensitymatrix!(ρ1, hMF(hartreefock), ks, hartreefock.μ; multimode=multimode, T=T, format=:dense) # get new meanfield and return the groundstate energy (density matrix was written to ρ1)
        # hartreefock(ρ0) # should actually not be needed here
        
        callback(ρ1)

        # if checkpoint != ""
        #     verbose ? @info("Saving intermediate mean field...") : nothing
        #     JLD.save(checkpoint, "mf", DenseHops(ρ1))
        # end

        ϵ0 + hartreefock.ϵMF # proper ground state energy
    end

    # Compute the ground state energy for the mean-field fixed point
    ϵ_GS, residual, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, verbose=verbose, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end
    hartreefock(ρ1) # update meanfield (h is updated in-place)


    dense(ρ1), ϵ_GS, hartreefock, converged, residual
end