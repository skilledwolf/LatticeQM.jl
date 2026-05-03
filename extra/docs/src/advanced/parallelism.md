# Parallelism

LatticeQM's parallel routines all go through one primitive — `Parallel.kspace_foreach!` for plain k-loops, `Parallel.kspace_reduce!` for additive accumulators (DOS, density, density-matrix, optical conductivity, Green's function). This page covers what knobs you have and when to turn them.

## Executor types

`Parallel` ships three executors:

| Executor | Activated by | When to use |
|----------|--------------|-------------|
| `SerialExec()` | always available | small problems, debugging, baseline measurements |
| `ThreadedExec(n; schedule=:dynamic)` | `julia -t N` or `JULIA_NUM_THREADS=N` | dense problems on a single machine; sparse moiré once you've moved off Arpack |
| `DistributedExec(n)` | `julia -p N` or `addprocs(N)` | multi-machine, very large k-grids, or sparse problems where threading isn't enough |

Pick one of three ways to specify the executor on any user-facing routine:

```julia
# Symbol shorthand — the historical multimode kwarg
bandmatrix(H, ks; multimode=:auto)            # default; picks the best available
bandmatrix(H, ks; multimode=:serial)
bandmatrix(H, ks; multimode=:threaded)        # alias :multithreaded
bandmatrix(H, ks; multimode=:distributed)     # gracefully degrades if no workers

# Or pass an Executor directly for fine control
bandmatrix(H, ks; executor=Parallel.ThreadedExec(8; schedule=:static))
bandmatrix(H, ks; executor=Parallel.DistributedExec(4))
```

`:auto` picks `DistributedExec` if `nworkers() > 1`, else `ThreadedExec` if `Threads.nthreads() > 1`, else `SerialExec`. `:distributed` falls back to threaded or serial when no workers are configured (preserves the historical `multimode=:distributed` semantics that downstream callers rely on).

## BLAS thread pinning

Each Julia worker (and each thread, via libblastrampoline) opens its own BLAS thread pool. On a 14-core box, `addprocs(4)` with default BLAS spawns roughly 4 × 14 = 56 contending BLAS threads — measured slowdowns of **100×** vs the same problem run with BLAS pinned to 1.

The first call to any non-serial-executor routine pins BLAS to 1 thread per worker via `Parallel.configure_blas!`. This is automatic; you don't need to call it yourself, but you can:

```julia
using Distributed
addprocs(4)
@everywhere using LatticeQM
Parallel.configure_blas!(Parallel.DistributedExec())   # pins BLAS=1 on every worker
```

If you *want* BLAS threading (e.g., one big serial dense diagonalisation per worker), set `BLAS.set_num_threads(...)` after the executor is activated.

## Building custom k-loops

If you need a parallelism pattern that's not already wrapped (e.g., a custom observable that doesn't fit `dos`/`density`/`Green`/`opticalconductivity`), use the primitives directly.

The examples below run as part of the doc build (a small graphene problem, no parallelism actually engaged — but the API exercise is real and they break the build if a future change drifts).

```@example parallelism
using LatticeQM
import LatticeQM.Parallel
import LatticeQM.Spectrum
import LatticeQM.Eigen

# Tiny test problem for the worked examples below.
lat = Geometries.honeycomb()
H   = Hops(); Operators.nearestneighbor!(H, lat)
ks  = LatticeQM.Structure.points(kpath(lat; num_points=16))
nothing # hide
```

### `kspace_foreach!` — plain iteration

```@example parallelism
exec = Parallel.to_executor(:serial)   # set to :auto in production code
Parallel.configure_blas!(exec; verbose=false)

# Output array (per-k results live here).
sum_of_negatives = zeros(Float64, size(ks, 2))

Parallel.kspace_foreach!(ks, exec;
    scratch_factory = () -> (Hcache = Spectrum.bloch_buffer(H, ks),),
) do scratch, j, k
    Hk = Spectrum.bloch!(scratch.Hcache, H, k)   # zero-alloc for AbstractHops + dense format
    ϵs, _ = Eigen.geteigen!(Hk)
    sum_of_negatives[j] = sum(ϵ for ϵ in real.(ϵs) if ϵ < 0)
end
sum_of_negatives[1:3]
```

`scratch_factory` is called **once per chunk** (one per thread / one per worker), not per k. Reuse heavy buffers (Hcache, U·D scratch, anything proportional to N²) here. Per-k allocations defeat the whole point of having parallelism on dense moiré.

### `kspace_reduce!` — additive accumulation

For observables where each k contributes additively to a global accumulator:

```@example parallelism
frequencies = LinRange(-3.5, 3.5, 200)
DOS = zeros(Float64, length(frequencies))

Parallel.kspace_reduce!(DOS, ks, exec;
    scratch_factory = () -> (Hcache = Spectrum.bloch_buffer(H, ks),),
) do local_dos, scratch, j, k
    Hk = Spectrum.bloch!(scratch.Hcache, H, k)
    ϵs = Eigen.geteigvals!(Hk)
    Spectrum.dos!(local_dos, real.(ϵs), frequencies; broadening=0.05)
end
DOS ./= size(ks, 2)
sum(DOS)  # ≈ 2 (orbitals/cell) × π / Γ-broadening — order-of-magnitude check
```

Each task accumulates into its own `local_dos` (lock-free), then `kspace_reduce!` folds the partials into `DOS` once at the end. The progress bar (when supplied via `progress=ProgressMeter.Progress(...)`) is funneled through a `Channel` (or `RemoteChannel` for distributed) so the bar isn't touched concurrently.

## Bloch-matrix helpers

Three helpers from `Spectrum` make in-place H(k) construction painless:

- `Spectrum.bloch_buffer(H, ks; format=:dense)` — preallocates a buffer for dense H(k); returns `nothing` for `format=:sparse` (sparse H(k) is built fresh per k).
- `Spectrum.bloch!(buf, H, k)` — fills `buf` with H(k). For `AbstractHops` + dense buffer it dispatches to the in-place `fouriersum!` (zero-alloc); for sparse / generic operators it falls back to `H(k)`.
- `Spectrum._build_H!(out, H, k)` — the lower-level in-place fill (extension point: define for your operator type to opt into the zero-alloc fast path).

For an `AbstractHops` Hamiltonian, the fast path saves ~1 GB of allocation churn / ~30% GC time on a TBG-N=5 bandmatrix call.

## Adding a custom output type

`kspace_reduce!` requires `output` to support `zero(output)` (for per-task accumulators) and `Parallel._accumulate!(out, partial)` (for the master-side fold). The default `_accumulate!` does broadcasted `out .+= partial`. If your output type doesn't broadcast cleanly, extend it:

```julia
# In your module:
import LatticeQM.Parallel
Base.zero(x::MyOutputType) = ...
Parallel._accumulate!(out::MyOutputType, partial::MyOutputType) = ...
```

`AbstractHops` does this — it opts out of broadcasting via `Base.broadcastable(H::Hops) = Ref(H)`, so `Operators` defines `_accumulate!` to walk δL keys and `Base.zero` to build a same-shape Hops with zero matrices.

## When threading vs distributed

Rough guidance based on measured benchmarks (Apple M4 Pro, 10 P-cores; results in `extra/benchmarks/results.csv`):

- **N ≲ 200 (dense)**: serial. Distributed loses to serial because of process startup; threading at best matches serial for trivial problems.
- **200 ≲ N ≲ 1000 (dense)**: threading. Distributed gives a tiny edge but at 4–5× the RSS.
- **Moiré (N ≳ 700, sparse)**: threading is the moiré default. KrylovKit Lanczos is thread-safe, so the per-process Hops copy that distributed pays is avoided. Use distributed only when threading isn't enough — typically meaning ≳ 1000 k-points or multi-node clusters.

The benchmark harness is reproducible: `bash extra/benchmarks/run_all.sh` regenerates the table on any machine.

## Multi-node / cluster jobs

When one machine isn't enough, the right shape is one process per node × N threads inside:

```julia
using Distributed
addprocs([("node1", 1), ("node2", 1), ...]; exeflags="--project=$(Base.active_project()) -t 16")
@everywhere using LatticeQM
@everywhere LinearAlgebra.BLAS.set_num_threads(1)

bandmatrix(H, ks; multimode=:distributed)
```

`SharedDenseHops` (the in-tree memory-sharing optimisation that fires for `DenseHops` under `multimode=:distributed`) only works on a single machine — for multi-node, every worker pays the H serialisation cost.

## Common pitfalls

- **Forgetting BLAS pinning**: Running `addprocs(4)` and `using LatticeQM` without touching BLAS gives the slowest possible configuration on a multi-core machine — measured at 100× slowdown vs single-threaded BLAS. The auto-pinning in `Parallel.configure_blas!` fires on first use of a parallel executor; if you have your own k-loops outside of `Parallel.kspace_*`, call it manually.
- **Sparse + threading without KrylovKit**: Arpack is not thread-safe. If you `include("eigen_sparse_arpack.jl")` (the legacy revert option in `Eigen.jl`), you must keep `multimode=:serial` or `:distributed` for sparse problems. The default KrylovKit-based solver is thread-safe.
- **Per-k allocations in custom k-loops**: If you write your own loop and allocate a 1000×1000 matrix per k, threading won't scale — Julia's GC is global and 8 threads producing garbage simultaneously contend on the same collector. Always preallocate via `scratch_factory`.
- **`multimode=:global`**: deprecated alias for `:auto`. The old `:global` resolved at module-load time; `:auto` re-evaluates each call so it sees workers added later.
