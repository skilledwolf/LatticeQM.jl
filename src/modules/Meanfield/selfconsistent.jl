# using JLD
import LatticeQM.Utils
import LatticeQM.Structure
import LatticeQM.Operators
import LatticeQM.Spectrum
import LatticeQM.TightBinding
using LatticeQM.Utils.Context

# Specialized cases
solvehartreefock(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v), filling, args...; kwargs...)
solvefock(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=true, fock=false), filling, args...; kwargs...)

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

using SharedArrays
import LatticeQM.Utils

solveselfconsistent!(::SerialContext, ρ0::T, ρ1::T, args...; kwargs...) where {T} = solveselfconsistent!(DummyContext(), Utils.dense(ρ0), Utils.dense(ρ1), args...; multimode=:serial, kwargs...)
function solveselfconsistent!(::DistributedContext, ρ0::T, ρ1::T, args...; kwargs...) where {T} 
    # workaround to fix memory allocation bug in julia 1.5 to 1.10
    # It is a somewhat hacky solution for now, but it reduces the memory allocation dramatically.
    # Before, every call to getdensitymatrix! would allocate a new shared array on disk for the density matrix,
    # bombarding both the file system and RAM, elluding the garbage collector.

    ρ0 = Utils.dense(ρ0)
    ρ1 = Utils.dense(ρ1)

    Ls = keys(ρ0)
    N = length(ρ0)
    Ns = size(ρ0)
    ρ0_backed = SharedArray(zeros(ComplexF64, Ns..., N))
    for (l_, L) in enumerate(Ls)
        ρ0_backed[:, :, l_] .= ρ0[L]
    end

    Ls = keys(ρ1)
    N = length(ρ1)
    Ns = size(ρ1)
    ρ1_backed = SharedArray(zeros(ComplexF64, Ns..., N))
    for (l_, L) in enumerate(Ls)
        ρ1_backed[:, :, l_] .= ρ1[L]
    end

    ρ0 = Hops(Dict(l => view(ρ0_backed, :, :, n_) for (n_, l) in enumerate(keys(ρ0))))
    ρ1 = Hops(Dict(l => view(ρ1_backed, :, :, n_) for (n_, l) in enumerate(keys(ρ1))))

    solveselfconsistent!(DummyContext(), ρ0, ρ1, args...; multimode=:distributed, kwargs...)
end

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
    # println("typeof hartreefock_object: ", typeof(hartreefock))
    # println("sparsity: ", sum(abs.(hartreefock.h[[0, 0]]) .> 1e-8) / length(hartreefock.h[[0, 0]]))

    update! = let hartreefock = hartreefock, ks = ks, filling = filling, T = T, multimode = multimode, verbose = verbose, callback = callback
        function f(ρ1, ρ0)
            hartreefock(ρ0) # update meanfield (h is updated in-place)
            println("sparsity: ", sum(abs.(hartreefock.hMF[[0, 0]]) .> 1e-9) / length(hartreefock.hMF[[0, 0]]))

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
        f
    end

    # function update!(ρ1, ρ0)
        
    #     hartreefock(ρ0) # update meanfield (h is updated in-place)

    #     verbose ? @info("Updating chemical potential for given filling...") : nothing
    #     hartreefock.μ = Spectrum.chemicalpotential(hMF(hartreefock), ks, filling; T=T, multimode=multimode)


    #     verbose ? @info("Updating the mean field density matrix...") : nothing
    #     ϵ0 = Operators.getdensitymatrix!(ρ1, hMF(hartreefock), ks, hartreefock.μ; multimode=multimode, T=T, format=:dense) # get new meanfield and return the groundstate energy (density matrix was written to ρ1)
    #     # hartreefock(ρ0) # should actually not be needed here
        
    #     callback(ρ1)

    #     # if checkpoint != ""
    #     #     verbose ? @info("Saving intermediate mean field...") : nothing
    #     #     JLD.save(checkpoint, "mf", DenseHops(ρ1))
    #     # end

    #     ϵ0 + hartreefock.ϵMF # proper ground state energy
    # end

    # Compute the ground state energy for the mean-field fixed point
    ϵ_GS, residual, converged = fixedpoint!(update!, ρ1, ρ0; iterations=iterations, tol=tol, verbose=verbose, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end
    hartreefock(ρ1) # update meanfield (h is updated in-place)


    Utils.dense(ρ1), ϵ_GS, hartreefock, converged, residual
end