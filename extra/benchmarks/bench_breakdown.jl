# Break down per-k cost: H build vs eigen vs misc.

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
using LatticeQM
import LatticeQM.TightBinding: fouriersum!
import LatticeQM.Eigen
import LatticeQM.Spectrum: handleprojector

case = ARGS[1]
if case == "tbg_n5"
    lat = Geometries.honeycomb_twisted(5)
elseif case == "tbg_n3"
    lat = Geometries.honeycomb_twisted(3)
end
T = Operators.graphene(lat)
ks = LatticeQM.Structure.points(kpath(lat; num_points=32))
N = LatticeQM.Structure.Lattices.countorbitals(lat)
projectors = handleprojector()

# warm up everything
Hwarm = Matrix{ComplexF64}(undef, N, N)
fouriersum!(Hwarm, T, ks[:, 1])
Eigen.geteigen!(copy(Hwarm))

# Measure phases separately, serial
H_time = 0.0
H_alloc = 0
eig_time = 0.0
eig_alloc = 0
n_iter = 5
for _ in 1:n_iter
    Hcache = Matrix{ComplexF64}(undef, N, N)
    GC.gc(); GC.gc()
    t1 = @timed for j in 1:size(ks, 2); fouriersum!(Hcache, T, ks[:, j]); end
    global H_time += t1.time; global H_alloc += t1.bytes
    GC.gc(); GC.gc()
    Hcache_filled = copy(Hcache)
    fouriersum!(Hcache, T, ks[:, 1])
    t2 = @timed for j in 1:size(ks, 2)
        copyto!(Hcache_filled, Hcache)
        Eigen.geteigen!(Hcache_filled)
    end
    global eig_time += t2.time; global eig_alloc += t2.bytes
end

println("# $case N=$N, 32 k-pts, averaged over $n_iter runs, BLAS=1")
println("H build (fouriersum!):  ", round(1000 * H_time / n_iter / 32; digits=3), " ms/k, ",
        round(H_alloc / n_iter / 32 / 1024; digits=1), " KB/k allocs")
println("eigen!:                 ", round(1000 * eig_time / n_iter / 32; digits=3), " ms/k, ",
        round(eig_alloc / n_iter / 32 / 1024; digits=1), " KB/k allocs")
println("total per k (serial):   ", round(1000 * (H_time + eig_time) / n_iter / 32; digits=3), " ms")
