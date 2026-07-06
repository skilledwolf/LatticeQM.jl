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

# Chern numbers from `getcherns` on the Haldane phase: the lower band must
# have C = +1 and the upper band C = -1 at ϕ = π/2 (sign flips at ϕ = -π/2).
# This is the only test that exercises Spectrum.statesgrid + Spectrum.berry
# end-to-end; without it, the k-loop refactor in berry.jl is unverified.
@testset "Haldane: Chern numbers via berry()" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    addhaldane_fast!(H, lat, x -> 0.1; ϕ=π/2)

    cherns = LatticeQM.Spectrum.getcherns(H, 30, 30)
    @test isapprox(cherns[1],  1.0; atol=0.05)
    @test isapprox(cherns[2], -1.0; atol=0.05)

    # Sign of ϕ flips the chirality.
    H2 = Hops()
    Operators.nearestneighbor!(H2, lat, -1.0)
    addhaldane_fast!(H2, lat, x -> 0.1; ϕ=-π/2)
    cherns2 = LatticeQM.Spectrum.getcherns(H2, 30, 30)
    @test isapprox(cherns2[1], -1.0; atol=0.05)
    @test isapprox(cherns2[2],  1.0; atol=0.05)
end

# gapcherns(): the same Haldane phase has a single gap at half filling. Its
# *cumulative* (non-abelian) Chern number — the manifold of the lower band —
# equals the lower-band Chern: +1 at ϕ=π/2, −1 at ϕ=−π/2. Exercises gap
# detection + Spectrum.berry over an occupied manifold end-to-end.
@testset "Haldane: gap Chern via gapcherns()" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    addhaldane_fast!(H, lat, x -> 0.1; ϕ=π/2)

    gaps = LatticeQM.Spectrum.gapcherns(H, 30; gaptol=1e-3)
    @test length(gaps) == 1
    @test gaps[1].n == 1
    @test gaps[1].chern == 1
    @test gaps[1].elo < 0 < gaps[1].ehi          # gap straddles E=0 (half filling)

    H2 = Hops()
    Operators.nearestneighbor!(H2, lat, -1.0)
    addhaldane_fast!(H2, lat, x -> 0.1; ϕ=-π/2)
    gaps2 = LatticeQM.Spectrum.gapcherns(H2, 30; gaptol=1e-3)
    @test length(gaps2) == 1
    @test gaps2[1].chern == -1
end

# Hofstadter Diophantine convention: with LatticeQM's peierlsoutplane flux
# sign and Spectrum.berry orientation, the exact gap Chern numbers obey
#   r = q·s − p·C   (ν = s − C·Φ),
# NOT r = q·s + p·C. This pins the relative sign of the Peierls and Berry
# conventions — a sign flip in either one breaks this test. The Haldane mass
# is essential: for plain NN honeycomb the π-flux (Φ=1/2) Dirac touchings
# masquerade as finite gaps on a coarse k-grid and violate the constraint.
@testset "Haldane: Hofstadter gap Cherns obey r = qs − pC" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    addhaldane_fast!(H, lat, x -> 0.1; ϕ=π/2)

    wd = Operators.hofstadter_cherns(H, lat, 5, 15; gaptol=1e-3)
    robust = findall(wd.width .>= 0.1)     # skip sub-grid pseudo-gaps
    @test length(robust) >= 40
    for i in robust
        p = numerator(wd.flux[i]); q = denominator(wd.flux[i])
        r = round(Int, wd.filling[i] * q)
        C = wd.chern[i]
        @test mod(r + p * C, q) == 0                       # r = q·s − p·C
        # the exact C must be a diophantine_cherns candidate whenever its
        # Středa intercept s falls inside the heuristic 0 ≤ s ≤ Nw window
        s = (r + p * C) ÷ q
        if 0 <= s <= 2
            @test C in last.(Operators.diophantine_cherns(p, q, r; Nw=2))
        end
    end

    # energies-only labeller: where its Diophantine label is unambiguous
    # (nsol == 1) it must equal the exact Berry-flux Chern, sign included.
    lab = Operators.hofstadter_gaplabels(H, lat, 5, 15; gaptol=1e-3, Nw=2)
    exact = Dict((wd.flux[i], round(Int, wd.filling[i] * denominator(wd.flux[i]))) => wd.chern[i]
                 for i in robust)
    nchecked = 0
    for j in eachindex(lab.flux)
        lab.nsol[j] == 1 || continue
        key = (lab.flux[j], lab.r[j])
        haskey(exact, key) || continue
        @test lab.chern[j] == exact[key]
        nchecked += 1
    end
    @test nchecked >= 15
end

# Known failure modes of the cheap Diophantine window, pinned as behavior.
# For plain NN honeycomb the gap at Φ = 4/5 with r = 3 has the exact label
# (s, C) = (3, 3) (r = q·s − p·C: 3 = 5·3 − 4·3): its Středa intercept s = 3
# lies OUTSIDE the naive 0 ≤ s ≤ Nw = 2 window, so smargin = 0 returns no
# admissible branch at all — the true branch only appears with smargin ≥ 1,
# and even then it is not the smallest-|C| candidate. This is why Level-1
# labels must be anchored (see hofstadter_gaplabels docstring).
@testset "diophantine_cherns: honeycomb out-of-window branch" begin
    @test isempty(Operators.diophantine_cherns(4, 5, 3; Nw=2))
    cands = Operators.diophantine_cherns(4, 5, 3; Nw=2, smargin=1)
    @test (3, 3) in cands
    @test first(cands) != (3, 3)          # smallest-|C| heuristic still wrong

    # all returned branches must actually solve r = q·s − p·C
    for (s, C) in cands
        @test 3 == 5 * s - 4 * C
    end

    # gaplabels must still REPORT such gaps (nsol == 0), not drop them
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    lab = Operators.hofstadter_gaplabels(H, lat, 5, 15; gaptol=1e-3, Nw=2)
    idx = findall(i -> lab.flux[i] == 4 // 5 && lab.r[i] == 3, eachindex(lab.flux))
    @test length(idx) == 1
    @test lab.nsol[idx[1]] == 0
end
