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
solvehartreefock(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent(ρ_init, HartreeFock(h, v), filling, args...; kwargs...)
solvefock(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=true, fock=false), filling, args...; kwargs...)

# Interface to solveselfconsistent!(ρ0, ...)
"""
    solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling, ks; kwargs...)
    solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling; klin, kwargs...)

Non‑mutating convenience wrappers around [`solveselfconsistent!`] that copy the
initial density matrix `ρ0`, iterate the mean‑field functional `mf` (e.g.,
`HartreeFock`), and return the converged result together with energy and state.

The `filling` sets the target electronic filling (0–1 per spin). Supply either
an explicit k‑grid `ks` or a grid resolution via `klin` (uses `klin×klin`).

Common keywords: `iterations`, `tol`, `T` (temperature), `β` (mixing),
`multimode` (parallel), `log_callback`.
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
function solveselfconsistent!(ρ0::T0, ρ1::T0, mf::MeanfieldGenerator, filling::Float64, ks; multimode=:auto, kwargs...) where {T0}
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

sanitize!(X) = X # dummy function, supply dispatch for your type
function sanitize!(ρ::TightBinding.Hops)
    if !TightBinding.ishermitian(ρ)
        @info "Initial guess is not hermitian, symmetrizing it now."
        TightBinding.hermitianize!(ρ)
    end
    ρ
end

"""
    _scf_driver!(ρ0, ρ1, mf, update_ρ!; kwargs...)

Common SCF skeleton shared by `solveselfconsistent!` (k-space diagonalization)
and `solveselfconsistent_purification!` (real-space canonical purification).
The caller supplies a closure `update_ρ!(ρ1, mf) -> ϵ_kinetic` that overwrites
`ρ1` from the current mean-field operator and returns the kinetic energy
contribution at this iteration.

Keywords:

  * `convergenceerror::Bool=false` — throw if SCF did not converge.
  * `callback=(ρ1 -> nothing)` — invoked with the freshly computed `ρ1` at
    each iteration. Existing user hook.
  * `log_callback=nothing` — if non-`nothing`, called as
    `log_callback(iter, ϵ, residual)` from within `fixedpoint!`. Lets callers
    log SCF progress without enabling the progress bar.
  * `iterations::Int=500`, `tol::Real=1e-7`, `verbose::Bool=false` — passed
    through to `fixedpoint!`.
  * Any other kwargs are forwarded to `fixedpoint!` (e.g. `β`, `p_norm`,
    `relative`, `show_trace`).
"""
function _scf_driver!(ρ0, ρ1, hartreefock, update_ρ!::F;
    convergenceerror=false, callback=(x -> nothing), log_callback=nothing,
    iterations=500, tol=1e-7, verbose::Bool=false, kwargs...) where {F}

    sanitize!(ρ0)
    sanitize!(ρ1)
    @assert TightBinding.ishermitian(ρ0) "Initial guess for density matrix must be hermitian."

    function update!(ρ1, ρ0)
        verbose ? @info("Updating mean field operators...") : nothing
        hartreefock(ρ0)

        verbose ? @info("Updating the mean field density matrix...") : nothing
        ϵ0 = update_ρ!(ρ1, hartreefock)

        callback(ρ1)
        ϵ0 + hartreefock.ϵMF
    end

    ϵ_GS, residual, converged = fixedpoint!(update!, ρ1, ρ0;
        iterations=iterations, tol=tol, verbose=verbose,
        log_callback=log_callback, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end
    hartreefock(ρ1)

    ρ_out = Utils.densecopy(ρ1)
    sanitize!(ρ_out)

    ρ_out, ϵ_GS, hartreefock, converged, residual
end

"""
    solveselfconsistent!(ρ0, ρ1, mf, filling, ks; kwargs...)
    solveselfconsistent!(ρ0, mf, filling, ks; kwargs...)
    solveselfconsistent!(ρ0, ρ1, mf, filling; klin, kwargs...)
    solveselfconsistent!(ρ0, mf, filling; klin, kwargs...)

Search for a self-consistent mean-field solution at the given `filling`
(between 0 and 1). The `MeanfieldGenerator` `mf` (e.g. `HartreeFock`) maps a
density matrix `ρ` to the effective single-particle Hamiltonian. k-space is
discretized with the supplied points `ks` (or built from `klin`).

Returns `(ρ, ϵ_GS, mf, converged, residual)`:

  1. converged density matrix,
  2. mean-field ground-state energy,
  3. the updated `MeanfieldGenerator` (with chemical potential `mf.μ`),
  4. convergence flag,
  5. final residual.

Keywords (all optional): `iterations`, `tol`, `T` (electronic temperature),
`β` (linear mixing), `convergenceerror`, `multimode`, `verbose`, `hidebar`,
`callback`, `log_callback`. See [`_scf_driver!`](@ref) and
[`fixedpoint!`](@ref) for details.

`parallel=true` (via `multimode`) helps when per-k diagonalization dominates
(e.g. twisted bilayer graphene). For small problems, the serial path can be
faster due to communication overhead.
"""
function _solveselfconsistent_impl!(ρ0::T1, ρ1::T1, hartreefock::MeanfieldGenerator, filling::Float64, ks;
    multimode=:serial, T=0.0, hidebar::Bool=false, verbose::Bool=false, kwargs...) where {T1}

    update_ρ!(ρ1, hf) = begin
        verbose ? @info("Updating chemical potential for given filling...") : nothing
        hf.μ = Spectrum.chemicalpotential(hMF(hf), ks, filling; T=T, multimode=multimode, hidebar=hidebar)
        ϵ0 = Operators.getdensitymatrix!(ρ1, hMF(hf), ks, hf.μ; multimode=multimode, T=T, format=:dense, hidebar=hidebar)
        if HARTREEFOCK_DEBUG[]
            @assert TightBinding.ishermitian(ρ1) "SANITY CHECK: HERMITIAN?"
        end
        ϵ0
    end

    _scf_driver!(ρ0, ρ1, hartreefock, update_ρ!; verbose=verbose, kwargs...)
end
