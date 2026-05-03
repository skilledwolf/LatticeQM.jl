# BLAS thread sweep: keep mode and case fixed, vary BLAS threads.
# Run as: julia --project=. bench_blas.jl <case> <mode> <blas_n>

using Distributed

case = ARGS[1]
mode = ARGS[2]
blas_n = parse(Int, ARGS[3])
nw = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 4

if mode == "distributed"
    addprocs(nw; exeflags=["--project=$(Base.active_project())"])
    @everywhere using LatticeQM
    @everywhere using LinearAlgebra
    @eval @everywhere LinearAlgebra.BLAS.set_num_threads($blas_n)
end

using LatticeQM
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(blas_n)

if case == "tbg_n7"
    lat = Geometries.honeycomb_twisted(7); ks = kpath(lat; num_points=32)
    extra = (; format=:sparse, num_bands=10)
elseif case == "tbg_n11"
    lat = Geometries.honeycomb_twisted(11); ks = kpath(lat; num_points=16)
    extra = (; format=:sparse, num_bands=10)
elseif case == "tbg_n5_dense"
    lat = Geometries.honeycomb_twisted(5); ks = kpath(lat; num_points=32)
    extra = (;)
end

T = Operators.graphene(lat)
multimode = mode == "threaded" ? :multithreaded :
            mode == "distributed" ? :distributed : :serial

# warm-up
warm_ks = kpath(lat; num_points=4)
Spectrum.bandmatrix(T, warm_ks; multimode=multimode, hidebar=true, extra...)

GC.gc(); GC.gc()
t = @timed Spectrum.bandmatrix(T, ks; multimode=multimode, hidebar=true, extra...)

println("$case,$mode,blas=$blas_n,nworkers=$(nworkers()),threads=$(Threads.nthreads()),time_s=$(round(t.time, digits=3))")
