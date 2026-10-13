import LatticeQM.Utils
import LatticeQM.Structure
import LatticeQM.Operators
import LatticeQM.Spectrum
import LatticeQM.TightBinding
import LatticeQM.Parallel

using LinearAlgebra: mul!
using Printf: @sprintf, @printf

# Specialized cases
"""
    solvehartreefock(h, v, ρ_init, filling; kwargs...)

Convenience wrapper around `solveselfconsistent` that constructs a `HartreeFock`
functional from base Hamiltonian `h` and interaction kernel `v`. Returns the
converged mean‑field solution and metadata.
"""
function solvehartreefock(h, v, ρ_init, filling::Number, args...; exchange=nothing, kwargs...)
    mf = exchange === nothing ? HartreeFock(h, v) : HartreeFock(h, v; exchange=exchange)
    solveselfconsistent(ρ_init, mf, filling, args...; kwargs...)
end
solvefock(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=false, fock=true), filling, args...; kwargs...)
solvehartree(h::T0, v, ρ_init, filling::Number, args...; kwargs...) where {T0} = solveselfconsistent(ρ_init, HartreeFock(h, v; hartree=true, fock=false), filling, args...; kwargs...)

# Interface to solveselfconsistent!(ρ0, ...)
"""
    solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling, ks; kwargs...)
    solveselfconsistent(ρ0, mf::MeanfieldGenerator, filling; klin, kwargs...)

Non‑mutating convenience wrappers around `solveselfconsistent!` that copy the
initial density matrix `ρ0`, iterate the mean‑field functional `mf` (e.g.,
`HartreeFock`), and return the converged result together with energy and state.

The `filling` sets the target electronic filling (0–1 per spin). Supply either
an explicit k‑grid `ks` or a grid resolution via `klin` (uses `klin×klin`).

Common keywords: `iterations`, `tol`, `T` (temperature), `β` (mixing),
`multimode` (parallel), `log_callback`, `residual_norm`.

# Returned energy decomposition

The returned `mf` carries the converged decomposition (physical-sign
convention):

  * `mf.ϵkin`  — kinetic energy `Tr[h₀ · ρ]` (derived from the variational
    identity `ϵkin = ϵband - 2(ϵH + ϵF)`, exact at the SCF fixed point)
  * `mf.ϵH`    — Hartree energy `+½ nᵀ V₀ n` (positive for repulsive `V₀`)
  * `mf.ϵF`    — Fock energy `-½ Σ_L Σ_{ij} v_L[i,j] |ρ_L[i,j]|²` (negative
    for repulsive `v`)
  * `mf.ϵband` — band energy `Tr[H_MF · ρ]` = sum of occupied eigenvalues

The total HF energy (returned as the second tuple element) can be expressed
in three equivalent ways at the SCF fixed point:

  * `E_HF = ϵkin + ϵH + ϵF`         (additive physical decomposition)
  * `E_HF = ϵband - ϵH - ϵF`        (band − double-counting)
  * `E_HF = ½ (ϵkin + ϵband)`       (`½ Tr[ρ (h₀ + H_MF)]`)

When `log_callback` is supplied, the driver invokes it as
`log_callback(iter, ϵ, residual, info)` if the callback can take a 4th arg
(detected via `applicable`), with `info = (ϵkin, ϵband, ϵH, ϵF, μ)`. Otherwise
it falls back to the legacy 3-arg form `log_callback(iter, ϵ, residual)`.
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

    # Worker reductions return ordinary dense matrices; SharedArrays created on
    # a worker need not be mapped on the process collecting the partials.
    mode = exec isa Parallel.DistributedExec ? :distributed :
           exec isa Parallel.ThreadedExec ? :multithreaded : :serial
    _solveselfconsistent_impl!(Utils.dense(ρ0), Utils.dense(ρ1), mf, filling, ks;
                              multimode=mode, kwargs...)
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
    commutator_kspace_norm(H, ρ, ks; multimode=:auto, kweights=nothing) → Real

Frobenius norm of the commutator `[H, ρ]`, evaluated in k-space:

```math
‖[H, ρ]‖_F = \\sqrt{\\sum_k w_k\\,‖H(k)\\,ρ(k) - ρ(k)\\,H(k)‖_F^2}
```

`H` and `ρ` must share an orbital basis (same matrix dimension per `Hops`
block); their key sets may differ. A vanishing commutator is necessary but
not sufficient for self-consistency: the occupations must also match the
requested filling and temperature. SCF combines this diagnostic with the
density fixed-point defect.
"""
function commutator_kspace_norm(H, ρ, ks::AbstractMatrix{Float64};
                                kweights::Union{Nothing,AbstractVector{Float64}}=nothing,
                                multimode::Symbol=:auto,
                                executor::Union{Nothing,Parallel.Executor}=nothing)
    nk = size(ks, 2)
    weights = kweights === nothing ? fill(1.0/nk, nk) : kweights

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    Parallel.configure_blas!(exec; verbose=false)

    factory = () -> begin
        Hbuf = Spectrum.bloch_buffer(H, ks; format=:dense)
        ρbuf = Spectrum.bloch_buffer(ρ, ks; format=:dense)
        N = size(Hbuf, 1)
        Cbuf = Matrix{ComplexF64}(undef, N, N)
        (Hcache=Hbuf, ρcache=ρbuf, C=Cbuf)
    end

    # `kspace_reduce!` requires an additive output, but we only need the
    # scalar reduction. A 1-element zero vector satisfies the contract;
    # `scalar_init=0.0` collects the per-k norm² contributions across tasks
    # / workers and returns their sum (distributed-safe — `Ref` capture
    # would silently lose worker writes).
    norm_sq = Parallel.kspace_reduce!([0.0], ks, exec;
        scratch_factory = factory,
        scalar_init = 0.0) do _dummy, scratch, j, k
        Hk = Spectrum.bloch!(scratch.Hcache, H, k)
        ρk = Spectrum.bloch!(scratch.ρcache, ρ, k)
        # C = Hk * ρk - ρk * Hk, computed via two BLAS3 calls into the
        # same buffer.
        mul!(scratch.C, Hk, ρk)
        mul!(scratch.C, ρk, Hk, -1, 1)
        weights[j] * sum(abs2, scratch.C)
    end

    sqrt(real(norm_sq))
end

# Wraps a user-supplied 3-arg `log_callback` so the SCF driver can pass a
# 4th NamedTuple of energy-decomposition fields without breaking older
# callbacks. `applicable` returns true for any signature that accepts the
# extra argument (concrete or untyped slot).
function _wrap_log_callback(cb, info_provider)
    cb === nothing && return nothing
    function (iter, ϵ, residual)
        info = info_provider()
        if applicable(cb, iter, ϵ, residual, info)
            cb(iter, ϵ, residual, info)
        else
            cb(iter, ϵ, residual)
        end
    end
end

# ----------------------------------------------------------------------------
# Tabular SCF logger (installed by `_scf_driver!` when verbose=true and the
# user hasn't supplied their own `log_callback`).
# ----------------------------------------------------------------------------

const _SCF_LOG_HEADER =
    "  iter        E_HF              E_kin              E_H              E_F          residual"
const _SCF_LOG_RULE = "  " * "─"^88

function _scf_log_print_header(io::IO=stderr)
    println(io)
    println(io, _SCF_LOG_RULE)
    println(io, _SCF_LOG_HEADER)
    println(io, _SCF_LOG_RULE)
end

function _scf_log_print_row(io::IO, iter::Integer, ϵ::Real, info::NamedTuple, residual::Real)
    @printf(io, "  %4d   %16.10f   %16.10f   %14.8f   %14.8f   %.3e\n",
            iter, ϵ, info.ϵkin, info.ϵH, info.ϵF, residual)
end

function _scf_log_print_footer(io::IO, iters::Integer, residual::Real, converged::Bool)
    println(io, _SCF_LOG_RULE)
    if converged
        @printf(io, "  SCF converged after %d iterations  (residual %.3e)\n", iters, residual)
    else
        @printf(io, "  SCF did NOT converge in %d iterations  (residual %.3e)\n", iters, residual)
    end
    println(io)
end

# Default tabular printer. Returns `(callback, iter_ref)` — the closure
# matches the 4-arg `log_callback` shape, and `iter_ref[]` exposes the
# last iteration count so the SCF driver can read it for the footer.
function _make_default_logger(io::IO)
    iter_ref = Ref(0)
    cb = function (iter::Integer, ϵ::Real, residual::Real, info::NamedTuple)
        _scf_log_print_row(io, iter, ϵ, info, residual)
        iter_ref[] = iter
    end
    cb, iter_ref
end

"""
    _scf_driver!(ρ0, ρ1, mf, update_ρ!; kwargs...)

SCF driver for k-space diagonalization.
The caller supplies a closure `update_ρ!(ρ1, mf) -> ϵ_band` that overwrites
`ρ1` from the current mean-field operator and returns the band-energy
contribution `Tr[H_MF · ρ_{k+1}]` at this iteration. The driver caches that
on `mf.ϵband` and derives `mf.ϵkin = ϵband - 2(ϵH + ϵF)` (variational
identity), so callers read the full physical decomposition
`(ϵkin, ϵH, ϵF, ϵband)` from the returned mean-field generator.

Keywords:

  * `convergenceerror::Bool=false` — throw if SCF did not converge.
  * `callback=(ρ1 -> nothing)` — invoked with the freshly computed `ρ1` at
    each iteration. Existing user hook.
  * `log_callback=nothing` — if non-`nothing`, called every iteration. The
    driver invokes it as either `log_callback(iter, ϵ, residual, info)` (4
    args) or the legacy `log_callback(iter, ϵ, residual)` (3 args), picked
    via `applicable`. `info` is a NamedTuple `(ϵkin, ϵband, ϵH, ϵF, μ)`.
  * `residual_fn=nothing` — optional additional convergence metric
    `(ρ1, ρ0) -> Real`. Both this metric and the density defect must meet
    the tolerance.
  * `iterations::Int=500`, `tol::Real=1e-7`, `verbose::Bool=false` — passed
    through to `fixedpoint!`.
  * Any other kwargs are forwarded to `fixedpoint!` (e.g. `β`, `p_norm`,
    `relative`, `show_trace`).
"""
function _scf_driver!(ρ0, ρ1, hartreefock, update_ρ!::F;
    convergenceerror=false, callback=(x -> nothing), log_callback=nothing,
    residual_fn=nothing,
    iterations=500, tol=1e-7, verbose::Bool=false,
    log_io::IO=stderr,
    kwargs...) where {F}

    sanitize!(ρ0)
    sanitize!(ρ1)
    @assert TightBinding.ishermitian(ρ0) "Initial guess for density matrix must be hermitian."

    # Default per-iteration logger: only installed when `verbose=true` and
    # the user hasn't supplied their own `log_callback`. Prints a clean
    # tabular row per iteration. Header/footer are printed by the driver
    # below (around `fixedpoint!`); `iter_ref` carries the last iteration
    # count back to the footer.
    is_default_logger = verbose && log_callback === nothing
    iter_ref = Ref(0)
    if is_default_logger
        log_callback, iter_ref = _make_default_logger(log_io)
    end

    function update!(ρ1, ρ0)
        hartreefock(ρ0)
        ϵband = real(update_ρ!(ρ1, hartreefock))
        hartreefock.ϵband = ϵband
        # Variational identity: `Tr[hMF·ρ] = ⟨T⟩ + 2·DC`. Derive `ϵkin`
        # rather than running a second Bloch loop just for `Tr[h·ρ]`.
        # Exact at the SCF fixed point; during iteration, picks up the same
        # ρ_old/ρ_new mismatch the band energy itself carries. DC is
        # ϵH + ϵF for plain HF and ϵH + ϵF + ϵP for BdG — omitting ϵP
        # overstated superconducting condensation energies by ½Σv|ρΔ|².
        hartreefock.ϵkin = ϵband - 2 * doublecounting(hartreefock)

        callback(ρ1)
        # E_GS = ϵband - DC (double-counting form)
        #      = ϵkin + DC  (additive physical form, equivalent)
        ϵband - doublecounting(hartreefock)
    end

    info_provider = () -> (
        ϵkin  = hartreefock.ϵkin,
        ϵband = hartreefock.ϵband,
        ϵH    = hartreefock.ϵH,
        ϵF    = hartreefock.ϵF,
        μ     = hartreefock.μ,
    )
    wrapped_log = _wrap_log_callback(log_callback, info_provider)

    is_default_logger && _scf_log_print_header(log_io)

    # When we own the per-iteration log, silence `fixedpoint!`'s own
    # final-summary `@info` (it would print after our footer). Its
    # non-convergence `@warn` still fires — we want that.
    fp_verbose = is_default_logger ? false : verbose
    ϵ_GS, residual, converged = fixedpoint!(update!, ρ1, ρ0;
        iterations=iterations, tol=tol, verbose=fp_verbose,
        log_callback=wrapped_log, residual_fn=residual_fn, kwargs...)

    if convergenceerror && !converged
        error("Convergence error.")
    end
    # Refresh `hMF, ϵH, ϵF` for the actually-converged ρ. `ϵband` still
    # reflects the last iteration's Bloch loop and we don't redo it (would
    # be a full extra `getdensitymatrix!`); re-derive `ϵkin` so the
    # variational identity `ϵkin = ϵband - 2(ϵH + ϵF)` holds exactly with
    # the current `ϵH, ϵF` values. At convergence ρ_{n-1} = ρ_n so the
    # returned `ϵ_GS` differs from `ϵband - ϵH - ϵF` only at SCF tolerance —
    # we replace it with the consistent value so the returned tuple matches
    # `mf` byte-for-byte.
    hartreefock(ρ1)
    hartreefock.ϵkin = hartreefock.ϵband - 2 * doublecounting(hartreefock)
    ϵ_GS = hartreefock.ϵband - doublecounting(hartreefock)

    is_default_logger && _scf_log_print_footer(log_io, iter_ref[], residual, converged)

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
    multimode=:serial, T=0.0, hidebar::Bool=true, verbose::Bool=false,
    residual_norm::Symbol=:density, kwargs...) where {T1}

    _ensure_density_support!(ρ0, hartreefock)
    _ensure_density_support!(ρ1, hartreefock)

    # `hidebar` defaults to `true` here (vs `false` in standalone calls):
    # the inner progress bars from `chemicalpotential` and
    # `getdensitymatrix!` fire once per SCF iteration, which makes the
    # screen unreadable. Set `hidebar=false` to re-enable them.
    update_ρ!(ρ1, hf) = begin
        hf.μ = Spectrum.chemicalpotential(hMF(hf), ks, filling; T=T, multimode=multimode, hidebar=hidebar)
        ϵ0 = Operators.getdensitymatrix!(ρ1, hMF(hf), ks, hf.μ; multimode=multimode, T=T, format=:dense, hidebar=hidebar)
        if HARTREEFOCK_DEBUG[]
            @assert TightBinding.ishermitian(ρ1) "SANITY CHECK: HERMITIAN?"
        end
        ϵ0
    end

    # The density defect is always checked before mixing. A commutator is
    # an additional condition: it alone cannot determine occupations. Its
    # Bloch reconstruction also requires enough density keys (enrichkeys!).
    residual_fn = if residual_norm === :density
        nothing
    elseif residual_norm === :commutator
        (ρ1_, ρ0_) -> commutator_kspace_norm(hMF(hartreefock), ρ0_, ks;
                                              multimode=multimode)
    else
        throw(ArgumentError("residual_norm must be :density or :commutator, got $(residual_norm)"))
    end

    _scf_driver!(ρ0, ρ1, hartreefock, update_ρ!;
                 verbose=verbose, residual_fn=residual_fn, kwargs...)
end
