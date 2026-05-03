# Compare LinearAlgebra.eigen! (current) against LAPACK.syevr!/heevr! with
# preallocated workspace buffers, to gauge whether workspace reuse closes the
# remaining threading gap on dense moiré.
#
# Hypothesis: per-k eigen allocates ~2.4 MB (the eigvec matrix) for N=364.
# 32 k × 8 threads × 2.4 MB = ~600 MB of garbage per bandmatrix call. If we
# can reuse the buffers, GC% drops further and threaded scaling improves.

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
using LatticeQM
import LatticeQM.TightBinding: fouriersum!

case = ARGS[1]
mode = ARGS[2]   # current | workspace

if case == "tbg_n5"
    lat = Geometries.honeycomb_twisted(5)
elseif case == "tbg_n3"
    lat = Geometries.honeycomb_twisted(3)
elseif case == "tbg_n7_dense"
    # force dense path even though autoconverter would pick sparse
    lat = Geometries.honeycomb_twisted(7)
end
T = Operators.graphene(lat)
ks = LatticeQM.Structure.points(kpath(lat; num_points=32))
N = LatticeQM.Structure.Lattices.countorbitals(lat)

# In-place Hermitian eigensolve using LAPACK.syevr! directly. Reuses the
# eigenvector buffer Z and the workspace W. For ComplexF64 input syevr! is
# the right call (despite the "sy" prefix it dispatches via libblastrampoline).
#
# Returns nothing — eigenvalues land in `vals`, eigenvectors in `Z`.
function eigen_inplace!(A::Matrix{ComplexF64}, vals::Vector{Float64}, Z::Matrix{ComplexF64})
    # syevr! signature on Complex: jobz='V' (compute eigvecs), range='A' (all),
    # uplo='U'. The work arrays are managed internally by LAPACK; what we
    # pass externally is just A (overwritten) and Z (output eigvecs).
    LAPACK.syevr!('V', 'A', 'U', A, 0.0, 0.0, 0, 0, -1.0)
end

# Warm-up: build one Hcache and invoke each path
Hwarm = Matrix{ComplexF64}(undef, N, N)
fouriersum!(Hwarm, T, ks[:, 1])

if mode == "current"
    # Current production code path
    function bench_current(N, ks, T)
        Hcache = Matrix{ComplexF64}(undef, N, N)
        for j in 1:size(ks, 2)
            fouriersum!(Hcache, T, ks[:, j])
            ϵs, U = eigen!(Hermitian(Hcache))   # what Julia stdlib does
        end
    end
    bench_current(N, ks, T)
    GC.gc(); GC.gc()
    t = @timed bench_current(N, ks, T)
    println("$case,current,time_s=$(round(t.time, digits=4)),alloc_mb=$(round(t.bytes/1024^2, digits=1)),gc_pct=$(round(100*t.gctime/t.time, digits=1))")

elseif mode == "workspace"
    # Hypothetical: reuse eigvec buffer
    function bench_workspace(N, ks, T)
        Hcache = Matrix{ComplexF64}(undef, N, N)
        # syevr! overwrites A (Hcache) in-place and returns (W::Vector{Float64}, Z::Matrix{ComplexF64})
        # The Z it returns is freshly allocated each call — we can't actually pre-pass it through
        # the Julia LAPACK wrapper. So this measures what's achievable WITHOUT touching Julia stdlib.
        for j in 1:size(ks, 2)
            fouriersum!(Hcache, T, ks[:, j])
            W, Z = LAPACK.syevr!('V', 'A', 'U', Hcache, 0.0, 0.0, 0, 0, -1.0)
        end
    end
    bench_workspace(N, ks, T)
    GC.gc(); GC.gc()
    t = @timed bench_workspace(N, ks, T)
    println("$case,workspace,time_s=$(round(t.time, digits=4)),alloc_mb=$(round(t.bytes/1024^2, digits=1)),gc_pct=$(round(100*t.gctime/t.time, digits=1))")

elseif mode == "current_threaded"
    function bench_t(N, ks, T)
        nt = Threads.nthreads()
        chunks = [1+(i-1)*cld(size(ks,2), nt):min(i*cld(size(ks,2), nt), size(ks,2)) for i in 1:nt]
        Threads.@threads for ic in 1:length(chunks)
            Hcache = Matrix{ComplexF64}(undef, N, N)
            for j in chunks[ic]
                fouriersum!(Hcache, T, ks[:, j])
                eigen!(Hermitian(Hcache))
            end
        end
    end
    bench_t(N, ks, T)
    GC.gc(); GC.gc()
    t = @timed bench_t(N, ks, T)
    println("$case,current_threaded,threads=$(Threads.nthreads()),time_s=$(round(t.time, digits=4)),alloc_mb=$(round(t.bytes/1024^2, digits=1)),gc_pct=$(round(100*t.gctime/t.time, digits=1))")

elseif mode == "workspace_threaded"
    function bench_tw(N, ks, T)
        nt = Threads.nthreads()
        chunks = [1+(i-1)*cld(size(ks,2), nt):min(i*cld(size(ks,2), nt), size(ks,2)) for i in 1:nt]
        Threads.@threads for ic in 1:length(chunks)
            Hcache = Matrix{ComplexF64}(undef, N, N)
            for j in chunks[ic]
                fouriersum!(Hcache, T, ks[:, j])
                LAPACK.syevr!('V', 'A', 'U', Hcache, 0.0, 0.0, 0, 0, -1.0)
            end
        end
    end
    bench_tw(N, ks, T)
    GC.gc(); GC.gc()
    t = @timed bench_tw(N, ks, T)
    println("$case,workspace_threaded,threads=$(Threads.nthreads()),time_s=$(round(t.time, digits=4)),alloc_mb=$(round(t.bytes/1024^2, digits=1)),gc_pct=$(round(100*t.gctime/t.time, digits=1))")
end
