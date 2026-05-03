using Test
using LatticeQM
using Distributed
import LinearAlgebra
import SparseArrays

# Parallel module: executor selection and the kspace_foreach! primitive.
@testset "Parallel: executor selection" begin
    @test Parallel.to_executor(:serial) isa Parallel.SerialExec
    # threaded with 1 thread degrades to serial (matches multimode=:multithreaded behaviour)
    if Threads.nthreads() > 1
        @test Parallel.to_executor(:threaded) isa Parallel.ThreadedExec
        @test Parallel.to_executor(:multithreaded) isa Parallel.ThreadedExec
    else
        @test Parallel.to_executor(:threaded) isa Parallel.SerialExec
    end
    # distributed degrades when no workers exist
    if nworkers() > 1
        @test Parallel.to_executor(:distributed) isa Parallel.DistributedExec
    else
        @test !(Parallel.to_executor(:distributed) isa Parallel.DistributedExec)
    end
    # auto picks the best available
    e = Parallel.to_executor(:auto)
    @test e isa Parallel.Executor
    @test Parallel.isavailable(e)
end

@testset "Parallel: kspace_foreach! correctness" begin
    ks = reshape(collect(0.0:0.1:1.0), 1, 11)  # 11 k-points
    # Per-thread/worker scratch should accumulate independently and produce the
    # same final answer as serial — body writes into a shared output array.
    for exec in (Parallel.SerialExec(),
                 Threads.nthreads() > 1 ? Parallel.ThreadedExec() : Parallel.SerialExec())
        results = zeros(Float64, size(ks, 2))
        Parallel.kspace_foreach!(ks, exec; scratch_factory = () -> nothing) do _, j, k
            results[j] = sin(k[1])  # idempotent per-j write — no race
        end
        @test results ≈ sin.(vec(ks))
    end
end

# Regression: `kspace_reduce!` on `ThreadedExec(>1)` was producing 1–10%
# non-deterministic drift on accumulators (filling/density/energy) before
# the closure-specialisation fix (`where {F, G, P}` on the threaded helper).
# This test reproduces the original drift signature on a small graphene
# density-matrix sweep and asserts the threaded result matches serial to
# floating-point precision across multiple runs.
@testset "Parallel: kspace_reduce! threaded matches serial" begin
    if Threads.nthreads() > 1
        lat = Geometries.honeycomb()
        h = Hops(); Operators.nearestneighbor!(h, lat, -1.0)
        ks = LatticeQM.Structure.points(LatticeQM.Structure.regulargrid(nk=12^2))

        ρS = Hops(); for k in keys(h); ρS[k] = zeros(ComplexF64, size(h[k])); end
        ϵS = Operators.getdensitymatrix!(ρS, h, ks, fill(1/size(ks,2), size(ks,2)), 0.0;
                                           T=0.05, hidebar=true, multimode=:serial)
        fS = real(LinearAlgebra.tr(ρS[[0, 0]])) / 2

        # 3 trials — drift was non-deterministic, would flake even with one
        # threaded run, but multiple trials make the test more sensitive.
        for _ in 1:3
            ρT = Hops(); for k in keys(h); ρT[k] = zeros(ComplexF64, size(h[k])); end
            ϵT = Operators.getdensitymatrix!(ρT, h, ks, fill(1/size(ks,2), size(ks,2)), 0.0;
                                               T=0.05, hidebar=true, multimode=:multithreaded)
            fT = real(LinearAlgebra.tr(ρT[[0, 0]])) / 2
            @test isapprox(fT, fS; atol=1e-10)
            @test isapprox(ϵT, ϵS; atol=1e-10)
        end
    end
end

# Regression: scalar accumulators captured by closure (Ref + SpinLock) used
# to be silently lost under `DistributedExec` because the closure was
# serialised per-worker. The fix is the `scalar_init` kwarg on
# `kspace_reduce!`, which collects body return values across tasks /
# workers. This test exercises serial / threaded / distributed and asserts
# the kinetic-energy scalar matches across all three.
@testset "Parallel: kspace_reduce! scalar_init across modes" begin
    lat = Geometries.honeycomb()
    h = Hops(); Operators.nearestneighbor!(h, lat, -1.0)
    ks = LatticeQM.Structure.points(LatticeQM.Structure.regulargrid(nk=8^2))
    nks = size(ks, 2)
    kweights = fill(1/nks, nks)

    ρS = Hops(); for k in keys(h); ρS[k] = zeros(ComplexF64, size(h[k])); end
    ϵS = Operators.getdensitymatrix!(ρS, h, ks, kweights, 0.0;
                                       T=0.05, hidebar=true, multimode=:serial)

    if Threads.nthreads() > 1
        ρT = Hops(); for k in keys(h); ρT[k] = zeros(ComplexF64, size(h[k])); end
        ϵT = Operators.getdensitymatrix!(ρT, h, ks, kweights, 0.0;
                                           T=0.05, hidebar=true, multimode=:multithreaded)
        @test isapprox(ϵT, ϵS; atol=1e-10)
    end

    # Distributed: skip if no extra workers; nworkers()==1 means we're
    # running with the master only (which would degrade to serial anyway).
    if nworkers() > 1
        ρD = Hops(); for k in keys(h); ρD[k] = zeros(ComplexF64, size(h[k])); end
        ϵD = Operators.getdensitymatrix!(ρD, h, ks, kweights, 0.0;
                                           T=0.05, hidebar=true, multimode=:distributed)
        @test isapprox(ϵD, ϵS; atol=1e-10)
        @test ϵD != 0.0  # silent loss bug would have returned 0
    end
end

# Each `multimode` for `bandmatrix` must produce identical bands and projector
# expectation values, otherwise we have a parallel-mode regression.
@testset "Spectrum: bandmatrix multimode equivalence" begin
    lat = Geometries.honeycomb()
    T = Hops()
    Operators.nearestneighbor!(T, lat)

    ks = kpath(lat; num_points=64)
    proj = Operators.valley(lat)

    bands_serial, obs_serial = Spectrum.bandmatrix(T, ks, proj; multimode=:serial, hidebar=true)

    @testset "multithreaded matches serial" begin
        if Threads.nthreads() > 1
            bands_mt, obs_mt = Spectrum.bandmatrix(T, ks, proj; multimode=:multithreaded, hidebar=true)
            @test bands_mt ≈ bands_serial
            @test obs_mt ≈ obs_serial
        else
            @test_skip Threads.nthreads() > 1
        end
    end

    @testset "distributed matches serial" begin
        if nworkers() > 1
            bands_d, obs_d = Spectrum.bandmatrix(T, ks, proj; multimode=:distributed, hidebar=true)
            @test bands_d ≈ bands_serial
            @test obs_d ≈ obs_serial
        else
            @test_skip nworkers() > 1
        end
    end
end

# Sparse + threaded was a forbidden combination under the Arpack default
# (Arpack is not thread-safe). With KrylovKit it works; this test pins the
# new behaviour.
@testset "Spectrum: sparse bandmatrix multimode equivalence" begin
    lat = Geometries.honeycomb_twisted(3)  # 148 orbitals, sparse-eligible
    T = Operators.graphene(lat)
    ks = kpath(lat; num_points=8)

    num_bands = 6
    bands_s, _ = Spectrum.bandmatrix(T, ks; multimode=:serial, format=:sparse,
                                     num_bands=num_bands, hidebar=true)

    if Threads.nthreads() > 1
        bands_t, _ = Spectrum.bandmatrix(T, ks; multimode=:multithreaded, format=:sparse,
                                         num_bands=num_bands, hidebar=true)
        # KrylovKit sorts by |λ-σ|, not by band index — sort each column for comparison.
        @test sort(bands_s; dims=1) ≈ sort(bands_t; dims=1) atol=1e-8
    else
        @test_skip Threads.nthreads() > 1
    end

    if nworkers() > 1
        bands_d, _ = Spectrum.bandmatrix(T, ks; multimode=:distributed, format=:sparse,
                                         num_bands=num_bands, hidebar=true)
        @test sort(bands_s; dims=1) ≈ sort(bands_d; dims=1) atol=1e-8
    else
        @test_skip nworkers() > 1
    end
end

# Regression for the Sn rewrite: it must build the same block-diagonal matrix
# as the previous `for i = 1:2:2N; mat[i:i+1, i:i+1] = σn` loop, i.e. spin is
# the FAST index, matching addspin / zeeman conventions.
@testset "Operators.Sn block-diagonal layout" begin
    lat = Geometries.honeycomb()
    N = LatticeQM.Structure.Lattices.countorbitals(lat)
    @assert N == 2

    σx = ComplexF64[0.0 1.0; 1.0 0.0]
    σy = ComplexF64[0.0 -1.0im; 1.0im 0.0]
    σz = ComplexF64[1.0 0.0; 0.0 -1.0]

    for (n, σ) in (([1.0, 0.0, 0.0], σx),
                   ([0.0, 1.0, 0.0], σy),
                   ([0.0, 0.0, 1.0], σz))
        sn = Operators.Sn(lat, n)
        ref = SparseArrays.spzeros(ComplexF64, 2N, 2N)
        for i in 1:2:2N
            ref[i:i+1, i:i+1] .= σ
        end
        @test sn ≈ ref
    end

    # Also verify the SparseMatrixCSC structure is preserved (the old @simd /
    # spzeros version returned a sparse matrix; a regression to dense would
    # blow up memory for large lattices).
    @test Operators.Sn(lat, [0.0, 0.0, 1.0]) isa SparseArrays.AbstractSparseMatrix
end
