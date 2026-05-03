using Test
using LatticeQM
using LatticeQM.Operators: addhaldane_naive!, addhaldane_fast!
using LinearAlgebra: eigvals

# Haldane model on honeycomb. Adds a complex next-nearest-neighbour hopping
# t₂ exp(±iϕ) with the chirality of the closing triangle, opening a mass gap
# at the Dirac points.
#
# Analytic full gap at K: 2·|m_H| = 6√3·|t₂·sin(ϕ)|.
# With t₂=0.1, ϕ=π/2 → gap = 6√3·0.1 ≈ 1.0392, so |E_K| = 3√3·t₂ ≈ 0.5196.
@testset "Haldane: honeycomb opens analytic gap at K" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    addhaldane_naive!(H, lat, x -> 0.1; ϕ=π/2)

    ev = sort(real.(eigvals(Matrix(H([1/3, -1/3])))))
    @test isapprox(ev[1], -3 * sqrt(3) * 0.1; atol=1e-12)
    @test isapprox(ev[2],  3 * sqrt(3) * 0.1; atol=1e-12)
end

# Regression: addhaldane_fast! must produce the same hopping matrix as the
# naive path. Pre-fix the fast path silently returned a zero (or doubled)
# Haldane term; this guards against regressions during cKDTree → KDTree
# refactors and any future spatial-lookup tweaks.
@testset "Haldane: fast path agrees with naive path" begin
    @testset "honeycomb (single cell)" begin
        lat = Geometries.honeycomb()
        Hn = Hops(); Operators.nearestneighbor!(Hn, lat, -1.0)
        Hf = Hops(); Operators.nearestneighbor!(Hf, lat, -1.0)
        addhaldane_naive!(Hn, lat, x -> 0.1; ϕ=π/2)
        addhaldane_fast!( Hf, lat, x -> 0.1; ϕ=π/2)

        # Same hopping vectors, same matrix entries.
        @test sort(collect(keys(Hn))) == sort(collect(keys(Hf)))
        for R in keys(Hn)
            @test Hn[R] ≈ Hf[R]
        end

        # And therefore the same spectrum at any k.
        for k in ([0.0, 0.0], [1/3, -1/3], [0.5, 0.5], [0.1, 0.2])
            evn = sort(real.(eigvals(Matrix(Hn(k)))))
            evf = sort(real.(eigvals(Matrix(Hf(k)))))
            @test maximum(abs, evn .- evf) < 1e-12
        end
    end

    @testset "honeycomb_AB bilayer" begin
        lat = Geometries.honeycomb_AB()
        Hn = Hops(); Operators.nearestneighbor!(Hn, lat, -1.0)
        Hf = Hops(); Operators.nearestneighbor!(Hf, lat, -1.0)
        addhaldane_naive!(Hn, lat, x -> 0.05; ϕ=π/2)
        addhaldane_fast!( Hf, lat, x -> 0.05; ϕ=π/2)

        evn = sort(real.(eigvals(Matrix(Hn([0.1, 0.2])))))
        evf = sort(real.(eigvals(Matrix(Hf([0.1, 0.2])))))
        @test maximum(abs, evn .- evf) < 1e-12
    end
end

# Pinned Haldane spectrum at K for honeycomb. Catches sign flips, missing
# √3 prefactor, and any change to the chirality convention.
@testset "Haldane: pinned spectrum at K" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    addhaldane_fast!(H, lat, x -> 0.1; ϕ=π/2)

    ev = sort(real.(eigvals(Matrix(H([1/3, -1/3])))))
    @test isapprox(ev[1], -0.5196152422706632; atol=1e-12)
    @test isapprox(ev[2],  0.5196152422706632; atol=1e-12)
end
