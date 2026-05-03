"""
    Parallel

Single primitive for k-space iteration across executors.

The package historically reimplemented `if multimode == :distributed && nprocs() > 1 ...`
inside every routine that loops over k-points. This module replaces that with
one parameterised loop:

    Parallel.kspace_foreach!(body!, ks, exec; scratch_factory, progress)

`body!(scratch, j, k)` runs once per k-index. `scratch` is a per-task buffer
(thread-local or worker-local), produced once per chunk by `scratch_factory()`,
not per k. This is the difference that lets dense bandmatrix on a 364-orbital
TBG cell drop from ~970 MB allocations to ~85 MB.

Supported executors:

- `SerialExec()`           — no parallelism.
- `ThreadedExec(n=Threads.nthreads(); schedule=:dynamic)` — spawns `n` tasks
  on the default thread pool, each consuming work via a shared `Channel`.
  Dynamic by default, so irregular per-k cost (e.g. sparse Krylov convergence
  varying near band touchings) is load-balanced.
- `DistributedExec(n=nworkers())` — uses `pmap` over chunks; one scratch per
  worker.

`configure_blas!(exec)` pins BLAS to one thread per worker for non-serial
executors. Skipping this is a 100× footgun on a multi-core box (each worker's
default BLAS pool oversubscribes the cores). Always call it once before
running the executor.
"""
module Parallel

using Distributed
using LinearAlgebra: BLAS
import ProgressMeter

export Executor, SerialExec, ThreadedExec, DistributedExec
export to_executor, isavailable, configure_blas!, kspace_foreach!, kspace_reduce!

abstract type Executor end

struct SerialExec <: Executor end

struct ThreadedExec <: Executor
    nthreads::Int
    schedule::Symbol  # :dynamic or :static
end
ThreadedExec(; schedule::Symbol=:dynamic) = ThreadedExec(Threads.nthreads(), schedule)
ThreadedExec(n::Int; schedule::Symbol=:dynamic) = ThreadedExec(n, schedule)

struct DistributedExec <: Executor
    nworkers::Int
end
DistributedExec() = DistributedExec(max(nworkers(), 1))

isavailable(::SerialExec) = true
isavailable(e::ThreadedExec) = e.nthreads >= 1 && Threads.nthreads() >= e.nthreads
isavailable(e::DistributedExec) = e.nworkers >= 1 && nworkers() >= e.nworkers

"""
    auto() -> Executor

Pick the most parallel executor that's actually available right now.
Distributed wins if workers exist, else threads, else serial.
"""
function auto()
    if nworkers() > 1
        return DistributedExec()
    elseif Threads.nthreads() > 1
        return ThreadedExec()
    else
        return SerialExec()
    end
end

"""
    to_executor(x) -> Executor

Coerce a user-facing value (`Executor`, `Symbol`, or `nothing`) to an
`Executor`. Symbols supported: `:serial`, `:multithreaded`/`:threaded`,
`:distributed`, `:auto`. Falls back to `SerialExec` with a warning for
unknown symbols (matches the previous `multimode` semantics).
"""
to_executor(e::Executor) = e
to_executor(::Nothing) = auto()
function to_executor(s::Symbol)
    if s === :serial
        SerialExec()
    elseif s === :multithreaded || s === :threaded
        Threads.nthreads() > 1 ? ThreadedExec() : SerialExec()
    elseif s === :distributed
        # Graceful degradation: fall back to threads or serial if no workers
        # are configured, matching the historical `multimode=:distributed`
        # semantics that downstream callers (e.g. chemicalpotential) rely on.
        if nworkers() > 1
            DistributedExec()
        elseif Threads.nthreads() > 1
            ThreadedExec()
        else
            SerialExec()
        end
    elseif s === :auto
        auto()
    else
        @warn "Unknown executor symbol $s; falling back to :serial"
        SerialExec()
    end
end

const _BLAS_CONFIGURED = Ref(false)

"""
    configure_blas!(exec; verbose=true) -> Int

For non-serial executors, set BLAS to 1 thread per worker (and on every
distributed worker). Idempotent: subsequent calls are a no-op.

Returns the BLAS thread count it ended up setting (or `BLAS.get_num_threads()`
if no change was needed).
"""
function configure_blas!(exec::Executor; verbose::Bool=true)
    if exec isa SerialExec
        return BLAS.get_num_threads()
    end
    _BLAS_CONFIGURED[] && return BLAS.get_num_threads()

    target = 1
    BLAS.set_num_threads(target)
    if exec isa DistributedExec && nworkers() > 1
        # Push BLAS=1 to every worker. `LinearAlgebra` may not be in `Main`
        # on the worker process, so import there first. Julia 1.12 world-age
        # semantics require `invokelatest` to use the freshly-imported binding
        # in the same expression — without it we get a noisy warning per
        # worker on every BLAS pinning.
        for w in workers()
            remotecall_wait(w) do
                Core.eval(Main, :(import LinearAlgebra))
                Base.invokelatest(() -> Main.LinearAlgebra.BLAS.set_num_threads(1))
            end
        end
    end
    _BLAS_CONFIGURED[] = true
    verbose && @info "Parallel: BLAS pinned to $target thread(s) per worker (executor=$(typeof(exec)))"
    return target
end

# ---------------------------------------------------------------------------
# kspace_foreach! - the only loop primitive
# ---------------------------------------------------------------------------

"""
    kspace_foreach!(body!, ks::AbstractMatrix, exec::Executor;
                    scratch_factory = () -> nothing,
                    progress = nothing)

Run `body!(scratch, j, @view ks[:, j])` for `j ∈ 1:size(ks, 2)` under the
given executor. `scratch` is task-local: `scratch_factory()` is called once
per chunk (per thread / per worker), never per k. `progress`, if given, is a
`ProgressMeter.Progress` — updates are funneled through a `Channel` so only
one task touches the bar.

Returns `nothing`. Side effects on whatever `body!` mutates (`bands`/`obs`
matrices, accumulated densities, etc.) are the caller's responsibility.
"""
function kspace_foreach!(body!, ks::AbstractMatrix, exec::Executor;
                          scratch_factory = () -> nothing,
                          progress::Union{Nothing,ProgressMeter.AbstractProgress}=nothing)
    nks = size(ks, 2)
    publish, finish = _progress_pump(progress, exec)
    try
        if exec isa SerialExec
            _kspace_serial!(body!, ks, nks, scratch_factory, publish)
        elseif exec isa ThreadedExec
            _kspace_threaded!(body!, ks, nks, scratch_factory, publish, exec)
        elseif exec isa DistributedExec
            _kspace_distributed!(body!, ks, nks, scratch_factory, publish, exec)
        else
            error("Unsupported executor $(typeof(exec))")
        end
    finally
        finish()
    end
    return nothing
end

function _kspace_serial!(body!, ks, nks, scratch_factory, publish)
    scratch = scratch_factory()
    @inbounds for j in 1:nks
        body!(scratch, j, view(ks, :, j))
        publish()
    end
end

function _kspace_threaded!(body!, ks, nks, scratch_factory, publish, exec::ThreadedExec)
    nt = max(1, exec.nthreads)

    if exec.schedule === :dynamic
        # Producer-consumer: a single shared work channel, each task pulls
        # the next k-index. Naturally load-balanced across irregular per-k
        # cost (sparse Krylov convergence, etc.).
        work_ch = Channel{Int}(nks)
        for j in 1:nks; put!(work_ch, j); end
        close(work_ch)
        @sync for _ in 1:nt
            Threads.@spawn begin
                scratch = scratch_factory()
                for j in work_ch
                    body!(scratch, j, view(ks, :, j))
                    publish()
                end
            end
        end
    else
        # Static: equal chunks, one scratch per task.
        chunks = _chunked_ranges(nks, nt)
        @sync for chunk in chunks
            Threads.@spawn begin
                scratch = scratch_factory()
                for j in chunk
                    body!(scratch, j, view(ks, :, j))
                    publish()
                end
            end
        end
    end
    return nothing
end

function _kspace_distributed!(body!, ks, nks, scratch_factory, publish, exec::DistributedExec)
    nw = max(1, exec.nworkers)
    chunks = _chunked_ranges(nks, nw)
    pmap(chunks) do chunk
        scratch = scratch_factory()
        for j in chunk
            body!(scratch, j, view(ks, :, j))
            publish()
        end
        return length(chunk)
    end
    return nothing
end

function _chunked_ranges(total::Int, nchunks::Int)
    sz = cld(total, nchunks)
    [(1 + (i-1)*sz):min(i*sz, total) for i in 1:nchunks if (i-1)*sz < total]
end

# How `kspace_reduce!` folds per-task accumulators into the master output.
# Default uses broadcasted `.+=` (works for arrays, complex vectors, etc).
# Modules can extend this for non-broadcastable types — TightBinding extends
# it for `AbstractHops` (which deliberately opts out of broadcasting via
# `Base.broadcastable(H::Hops) = Ref(H)`, so `.+=` would fail).
@inline _accumulate!(out, partial) = (out .+= partial; out)

# Progress pump abstraction used by both kspace_foreach! and kspace_reduce!
# under threaded / distributed executors. Returns `(publish, finish)`:
#   - `publish()` is called by every per-k body iteration; thread-safe
#   - `finish()` is called once after all bodies complete; drains the pump
# For `progress === nothing` both are no-ops. For SerialExec, `publish` calls
# `next!` directly (no Channel overhead) since there's no race.
function _progress_pump(progress, exec::Executor)
    if progress === nothing
        return (() -> nothing), (() -> nothing)
    end
    if exec isa SerialExec
        return (() -> ProgressMeter.next!(progress)), (() -> nothing)
    end
    ch = exec isa DistributedExec ? RemoteChannel(() -> Channel{Bool}(Inf), 1) :
                                    Channel{Bool}(Inf)
    pump = @async while take!(ch)
        ProgressMeter.next!(progress)
    end
    publish = () -> put!(ch, true)
    finish = () -> (put!(ch, false); wait(pump))
    return publish, finish
end

# ---------------------------------------------------------------------------
# kspace_reduce! - additive reduction over k-points
# ---------------------------------------------------------------------------

"""
    kspace_reduce!(body!, output, ks, exec; scratch_factory, progress=nothing)

Run `body!(local_out, scratch, j, k)` once per k-index. `local_out` is a
per-task accumulator (one per chunk on threaded/distributed; just `output`
itself on serial). After the loop, `local_out`s are summed into `output` via
`output .+= partial`.

This is the right pattern for DOS, density, density-matrix, optical
conductivity — anywhere each k contributes additively to a global
accumulator. It is lock-free in the hot loop: each task writes only to its
own `local_out`, and the master merges them once at the end.

`progress`, if given, is a `ProgressMeter.Progress` — updates are funneled
through a `Channel` so only one task touches the bar.

`output` must support `zero(output)` and broadcasted `.+=`.

Returns `output`.
"""
function kspace_reduce!(body!, output, ks::AbstractMatrix, exec::Executor;
                         scratch_factory = () -> nothing,
                         progress::Union{Nothing,ProgressMeter.AbstractProgress}=nothing)
    nks = size(ks, 2)
    publish, finish = _progress_pump(progress, exec)

    try
        if exec isa SerialExec
            scratch = scratch_factory()
            for j in 1:nks
                body!(output, scratch, j, view(ks, :, j))
                publish()
            end
        elseif exec isa ThreadedExec
            chunks = _chunked_ranges(nks, exec.nthreads)
            partials = [zero(output) for _ in 1:length(chunks)]
            @sync for (ic, chunk) in enumerate(chunks)
                Threads.@spawn begin
                    scratch = scratch_factory()
                    for j in chunk
                        body!(partials[ic], scratch, j, view(ks, :, j))
                        publish()
                    end
                end
            end
            for p in partials
                _accumulate!(output, p)
            end
        elseif exec isa DistributedExec
            chunks = _chunked_ranges(nks, exec.nworkers)
            partials = pmap(chunks) do chunk
                local_out = zero(output)
                scratch = scratch_factory()
                for j in chunk
                    body!(local_out, scratch, j, view(ks, :, j))
                    publish()
                end
                return local_out
            end
            for p in partials
                _accumulate!(output, p)
            end
        else
            error("Unsupported executor $(typeof(exec))")
        end
    finally
        finish()
    end
    return output
end

end # module Parallel
