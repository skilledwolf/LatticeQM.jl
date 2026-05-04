# using JLD
import LatticeQM.Utils
import LatticeQM.Structure
import LatticeQM.Operators
import LatticeQM.Spectrum
import LatticeQM.TightBinding
import LatticeQM.Parallel

# Specialized cases
"""
    solvehartreefock(h, v, ρ_init, filling; kwargs...)

Convenience wrapper around `solveselfconsistent` that constructs a `HartreeFock`
functional from base Hamiltonian `h` and interaction kernel `v`. Returns the
converged mean‑field solution and metadata.
"""
solvehartreefock(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v), filling, args...; kwargs...)
solvefock(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree(h::T, v, ρ_init, filling::Number, args...; kwargs...) where {T} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=true, fock=false), filling, args...; kwargs...)

# Interface to solveselfconsistent!(ρ0, ...)
"""
    solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling, ks; kwargs...)
    solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling; klin, kwargs...)

Non‑mutating convenience wrappers around [`solveselfconsistent!`] that copy the
initial density matrix `ρ0`, iterate the mean‑field functional `mf` (e.g.,
`HartreeFock`), and return the converged result together with energy and state.

The `filling` sets the target electronic filling (0–1 per spin). Supply either
an explicit k‑grid `ks` or a grid resolution via `klin` (uses `klin×klin`).

Common keywords: `iterations`, `tol`, `T`, `β` (mixing), `multimode` (parallel).
"""
solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling::Number, ks; kwargs...) = solveselfconsistent!(deepcopy(ρ0), mf, filling, ks; kwargs...)
solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling::Number; klin, kwargs...) = solveselfconsistent!(deepcopy(ρ0), mf, filling; klin=klin, kwargs...)

# Interface to solveselfconsistent!(ρ0, ρ1, ...)
function solveselfconsistent!(ρ0, mf::MeanfieldGenerator, filling::Float64, args...; kwargs...)
    sanitize!(ρ0)
    solveselfconsistent!(ρ0, deepcopy(ρ0), mf, filling, args...; kwargs...)
end

# Interface to translate klin to solveselfconsistent!(ρ0, ρ1, ..., ks, ...)
solveselfconsistent!(ρ0, ρ1, mf::MeanfieldGenerator, filling::Float64; klin::Int, kwargs...) = solveselfconsistent!(ρ0, ρ1, mf, filling, Structure.regulargrid(nk=klin^2); kwargs...)

# `multimode=:global` is preserved for backwards compatibility but now maps
# to `:auto` (the old `:global` resolved to `getautocontext()` *at module-load
# time*, which froze to whatever Julia was started with — a footgun if
# `addprocs` happened later). `:auto` re-evaluates each call.
function solveselfconsistent!(ρ0::T, ρ1::T, mf::MeanfieldGenerator, filling::Float64, ks; multimode=:auto, kwargs...) where {T}
    mm = (multimode === :global) ? :auto : multimode
    exec = Parallel.to_executor(mm)

    if exec isa Parallel.DistributedExec
        # Convert to SharedDenseHops so the in-tree
        # `Spectrum.sanatize_distributed_hamiltonian(::DenseHops)` hook can hand
        # workers a shared buffer instead of full per-worker copies. The
        # historical extra step of wrapping each matrix in a `view` (via
        # `gethopsview` → `SubarrayHops`) is no longer needed: the per-task
        # accumulator is allocated by `Parallel.kspace_reduce!` via
        # `Base.zero(::Hops)`, so workers don't write into the SharedMatrix
        # directly anymore.
        return _solveselfconsistent_impl!(TightBinding.shareddense(ρ0),
                                           TightBinding.shareddense(ρ1),
                                           mf, filling, ks;
                                           multimode=:distributed, kwargs...)
    elseif exec isa Parallel.ThreadedExec
        # Threading shares memory by reference; a plain dense copy is fine.
        # (Previously this path errored — the legacy Context dispatch had no
        # MultiThreadedContext implementation. With Parallel.kspace_* under
        # the hood, threading just works.)
        return _solveselfconsistent_impl!(Utils.dense(ρ0), Utils.dense(ρ1), mf, filling, ks; multimode=:multithreaded, kwargs...)
    else  # SerialExec
        return _solveselfconsistent_impl!(Utils.dense(ρ0), Utils.dense(ρ1), mf, filling, ks; multimode=:serial, kwargs...)
    end
end

using SharedArrays

sanitize!(X) = X # dummy function, supply dispatch for your type
function sanitize!(ρ::TightBinding.Hops)
    if !TightBinding.ishermitian(ρ)
        @info "Initial guess is not hermitian, symmetrizing it now."
        TightBinding.hermitianize!(ρ)
    end
    ρ
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
function _solveselfconsistent_impl!(ρ0::T1, ρ1::T1, hartreefock::MeanfieldGenerator, filling::Float64, ks;
    convergenceerror=false, multimode=:serial, checkpoint::String="", callback=(x -> nothing), hotstart=true, iterations=500, tol=1e-7, T=0.0, format=:dense, verbose::Bool=false,
    hidebar::Bool=false, kwargs...) where {T1}

    # if checkpoint != "" && isfile(checkpoint) && hotstart
    #     println("Loading checkpoint file as initial guess: $checkpoint")
    #     ρ0 = JLD.load(checkpoint, "mf")
    # end

    sanitize!(ρ0)
    sanitize!(ρ1)
    @assert TightBinding.ishermitian(ρ0) "Initial guess for density matrix must be hermitian."

    function update!(ρ1, ρ0)
        verbose ? @info("Updating mean field operators...") : nothing
        hartreefock(ρ0) # update meanfield (h is updated in-place)
        # println("sparsity: ", sum(abs.(hartreefock.hMF[[0, 0]]) .> 1e-9) / length(hartreefock.hMF[[0, 0]]))

        verbose ? @info("Updating chemical potential for given filling...") : nothing
        hartreefock.μ = Spectrum.chemicalpotential(hMF(hartreefock), ks, filling; T=T, multimode=multimode, hidebar=hidebar)

        verbose ? @info("Updating the mean field density matrix...") : nothing
        ϵ0 = Operators.getdensitymatrix!(ρ1, hMF(hartreefock), ks, hartreefock.μ; multimode=multimode, T=T, format=:dense, hidebar=hidebar) # get new meanfield and return the groundstate energy (density matrix was written to ρ1)

        @assert TightBinding.ishermitian(ρ1) "SANITY CHECK: HERMITIAN?"

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

    ρout = Utils.densecopy(ρ1)
    sanitize!(ρout)

    ρout, ϵ_GS, hartreefock, converged, residual
end 
