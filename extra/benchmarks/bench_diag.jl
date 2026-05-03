# Diagnostic: where does dense bandmatrix time and memory go?

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
using LatticeQM
import LatticeQM.TightBinding: fouriersum!
import LatticeQM.Eigen
import LatticeQM.Spectrum: bandmatrix_Hcache, bandmatrix_Ucache,
    bandmatrix_preallocate, handleprojector, insertbands_bandexpvals_k!

case = ARGS[1]
mode = ARGS[2]   # "current" | "inplace"

if case == "tbg_n3"
    lat = Geometries.honeycomb_twisted(3)
elseif case == "tbg_n5"
    lat = Geometries.honeycomb_twisted(5)
end
T = Operators.graphene(lat)
ks = LatticeQM.Structure.points(kpath(lat; num_points=32))
N = LatticeQM.Structure.Lattices.countorbitals(lat)
nhops = length(T.data)
println("# orbitals=$N, hop offsets=$nhops, k-points=$(size(ks, 2))")

projectors = handleprojector()
bands, obs = bandmatrix_preallocate(T, ks, projectors)

# warm-up
Hwarm = Matrix{ComplexF64}(undef, N, N)
fouriersum!(Hwarm, T, ks[:, 1])

if mode == "current"
    # exactly what bandmatrix_serial! does:
    GC.gc(); GC.gc()
    t = @timed begin
        Hcache = Matrix{ComplexF64}(undef, N, N)
        for j in 1:size(ks, 2)
            Hcache .= T(ks[:, j])  # T(k) allocates ~nhops intermediate matrices
            ϵs, U = Eigen.geteigen!(Hcache)
            @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], ϵs, U, ks[:, j], projectors)
        end
    end
    println("current,time_s=", round(t.time, digits=3),
            ",alloc_mb=", round(t.bytes/1024^2, digits=1),
            ",gctime_s=", round(t.gctime, digits=3),
            ",gc_pct=", round(100*t.gctime/t.time, digits=1))

elseif mode == "inplace"
    # what the rewrite would do:
    GC.gc(); GC.gc()
    t = @timed begin
        Hcache = Matrix{ComplexF64}(undef, N, N)
        for j in 1:size(ks, 2)
            fouriersum!(Hcache, T, ks[:, j])  # in-place: zero alloc per call
            ϵs, U = Eigen.geteigen!(Hcache)
            @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], ϵs, U, ks[:, j], projectors)
        end
    end
    println("inplace,time_s=", round(t.time, digits=3),
            ",alloc_mb=", round(t.bytes/1024^2, digits=1),
            ",gctime_s=", round(t.gctime, digits=3),
            ",gc_pct=", round(100*t.gctime/t.time, digits=1))

elseif mode == "inplace_threaded"
    nthr = Threads.nthreads()
    chunks = [1+(i-1)*cld(size(ks,2), nthr):min(i*cld(size(ks,2), nthr), size(ks,2)) for i in 1:nthr]
    GC.gc(); GC.gc()
    t = @timed begin
        Threads.@threads for ic in 1:length(chunks)
            chunk = chunks[ic]
            Hloc = Matrix{ComplexF64}(undef, N, N)
            for j in chunk
                fouriersum!(Hloc, T, ks[:, j])
                ϵs, U = Eigen.geteigen!(Hloc)
                @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], ϵs, U, ks[:, j], projectors)
            end
        end
    end
    println("inplace_threaded,threads=$nthr,time_s=", round(t.time, digits=3),
            ",alloc_mb=", round(t.bytes/1024^2, digits=1),
            ",gctime_s=", round(t.gctime, digits=3),
            ",gc_pct=", round(100*t.gctime/t.time, digits=1))

elseif mode == "current_threaded"
    nthr = Threads.nthreads()
    chunks = [1+(i-1)*cld(size(ks,2), nthr):min(i*cld(size(ks,2), nthr), size(ks,2)) for i in 1:nthr]
    GC.gc(); GC.gc()
    t = @timed begin
        Threads.@threads for ic in 1:length(chunks)
            chunk = chunks[ic]
            Hloc = Matrix{ComplexF64}(undef, N, N)
            for j in chunk
                Hloc .= T(ks[:, j])
                ϵs, U = Eigen.geteigen!(Hloc)
                @views insertbands_bandexpvals_k!(bands[:, j], obs[:, j, :], ϵs, U, ks[:, j], projectors)
            end
        end
    end
    println("current_threaded,threads=$nthr,time_s=", round(t.time, digits=3),
            ",alloc_mb=", round(t.bytes/1024^2, digits=1),
            ",gctime_s=", round(t.gctime, digits=3),
            ",gc_pct=", round(100*t.gctime/t.time, digits=1))
end
