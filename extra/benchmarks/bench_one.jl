# Single-shot parallelism benchmark.
#
# Run as:  julia --project=. extra/benchmarks/bench_one.jl <case> <mode> [extra]
#
# <case> ∈ {dense_small, tbg_n3, tbg_n5, tbg_n7, tbg_n11}
# <mode> ∈ {serial, threaded, distributed}
#
# Output: one CSV line on stdout: case,mode,nthreads,nworkers,blas,nks,N,time_s,rss_mb,allocs_mb
#
# Wall time covers ONE bandmatrix call after a warm-up call.
# RSS is Sys.maxrss() at the end of the run (peak resident set since process start).
# Allocations are counted by @timed on the timed call.

using Distributed

case = ARGS[1]
mode = ARGS[2]

if mode == "distributed"
    nw = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 4
    addprocs(nw; exeflags=["--project=$(Base.active_project())"])
    @everywhere using LatticeQM
    @everywhere using LinearAlgebra
    @everywhere LinearAlgebra.BLAS.set_num_threads(1)
end

using LatticeQM
using LinearAlgebra
import Statistics
LinearAlgebra.BLAS.set_num_threads(1)  # avoid oversubscription

build_problem = function(case)
    if case == "dense_small"
        lat = Geometries.honeycomb()
        T = Hops()
        Operators.nearestneighbor!(T, lat)
        ks = kpath(lat; num_points=128)
        return (lat, T, ks, (;), 2)
    end
    n = parse(Int, replace(case, "tbg_n" => ""))
    lat = Geometries.honeycomb_twisted(n)
    T = Operators.graphene(lat)  # builds NN + interlayer hops, auto dense/sparse
    norbital = LatticeQM.Structure.Lattices.countorbitals(lat)
    nks = norbital > 700 ? 16 : 32
    ks = kpath(lat; num_points=nks)
    # for sparse path we need num_bands; pick a few low-energy bands
    extra = norbital > 500 ? (; format=:sparse, num_bands=10) : (;)
    return (lat, T, ks, extra, norbital)
end

lat, T, ks, extra, norbital = build_problem(case)

multimode = mode == "threaded" ? :multithreaded :
            mode == "distributed" ? :distributed : :serial

# warm-up (precompile, cache, etc.)
warm_ks = kpath(lat; num_points=4)
Spectrum.bandmatrix(T, warm_ks; multimode=multimode, hidebar=true, extra...)

GC.gc(); GC.gc()
t = @timed Spectrum.bandmatrix(T, ks; multimode=multimode, hidebar=true, extra...)
elapsed = t.time
allocs_mb = t.bytes / 1024^2
rss_mb = Sys.maxrss() / 1024^2
nks = size(LatticeQM.Structure.points(ks), 2)

# child workers' RSS is not counted in master's maxrss; sum them in distributed mode
total_rss_mb = rss_mb
if mode == "distributed"
    worker_rss = sum(remotecall_fetch(() -> Sys.maxrss()/1024^2, w) for w in workers())
    total_rss_mb = rss_mb + worker_rss
end

println(join((case, mode, Threads.nthreads(), nworkers(), LinearAlgebra.BLAS.get_num_threads(),
              nks, norbital, round(elapsed, digits=3), round(total_rss_mb, digits=1),
              round(allocs_mb, digits=1)), ","))
