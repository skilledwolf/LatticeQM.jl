using Test
using LinearAlgebra
using LatticeQM
using LatticeQM.Spectrum: statesgrid1D, Wilson1D, wilsonweights, sectorwilson,
                          WilsonWeights1D, SectorWilson1D, statesgrid1Dslice

# Two-leg ladder of SSH chains, basis (A1, B1, A2, B2), rung tperp on A1-B2,
# bias ±V. Leg ℓ has intracell amplitude tp[ℓ] and intercell amplitude t = 1.
function ladder(tp1, tp2, tperp, V)
    h = DenseHops()
    h[[0]] = ComplexF64[V tp1 0 0; tp1 V 0 0; 0 0 -V tp2; 0 0 tp2 -V]
    h[[1]] = zeros(ComplexF64, 4, 4)
    h[[1]][1, 2] = 1.0         # A1(n) -> B1(n+1)
    h[[1]][3, 4] = 1.0         # A2(n) -> B2(n+1)
    h[[1]][1, 4] = tperp       # A1(n) -> B2(n+1) (Bernal rung)
    h[[-1]] = h[[1]]'
    h
end
τz = Diagonal([1.0, 1.0, -1.0, -1.0])

@testset "wilsonweights and sectorwilson: decoupled legs" begin
    # leg 1 topological (|tp| < t), leg 2 trivial, no rung, no bias
    h = ladder(0.4, 1.6, 0.0, 0.0)
    centres, w = WilsonWeights1D(h, 201, Matrix(τz), 1:2)
    @test sort(round.(mod.(centres .+ 0.25, 1) .- 0.25; digits=6)) == [0.0, 0.5]
    # the half-integer centre lives on leg 1 (τz = +1)
    i = argmax(abs.(mod.(centres .+ 0.25, 1) .- 0.25))
    @test isapprox(w[i], 1.0; atol=1e-6)
    p, gap = SectorWilson1D(h, 201, Matrix(τz), 1:2)
    @test isapprox(gap, 1.0; atol=1e-9)
    # sector + (second entry, positive eigenvalue) carries 1/2, sector - carries 0
    @test isapprox(abs(p[2]), 0.5; atol=1e-6)
    @test isapprox(p[1], 0.0; atol=1e-6)
    # with decoupled legs the two constructions agree
    @test isapprox(sum(w .* centres) - (0.5 * 1.0), 0.0; atol=1e-6) || isapprox(abs(sum(w .* centres)), 0.5; atol=1e-6)
end

@testset "sectorwilson: coupled legs, gauge invariance, rung gap" begin
    h = ladder(0.4, 1.6, 0.3, 0.15)
    p, gap = SectorWilson1D(h, 301, Matrix(τz), 1:2)
    @test 0 < gap <= 1
    # phases are mod 1 and change by less than 1e-3 when the grid is refined
    p2, _ = SectorWilson1D(h, 601, Matrix(τz), 1:2)
    @test maximum(abs.(mod.(p .- p2 .+ 0.5, 1) .- 0.5)) < 1e-3
    # the holonomy-first weights lie in [-1, 1] and the centres are mod 1
    centres, w = WilsonWeights1D(h, 301, Matrix(τz), 1:2)
    @test all(abs.(w) .<= 1 + 1e-12)
    @test all(-0.5 .<= centres .<= 0.5)
end

@testset "statesgrid1Dslice: 2D Hamiltonian at fixed momentum" begin
    lat = Geometries.honeycomb()
    H = Operators.graphene(lat)
    _, sg = statesgrid1Dslice(H, 101, [1], 1, 0.2)   # loop along k2 at k1 = 0.2
    @test size(sg) == (101, 2, 1)
    _, sg2 = statesgrid1Dslice(H, 101, [1], 2, 0.2)  # loop along k1 at k2 = 0.2
    @test size(sg2) == (101, 2, 1)
    # Zak phase of the lower graphene band along k2 at fixed k1: 0 or 1/2 mod 1
    # (chiral symmetry); away from the Dirac projections it is well defined.
    z = Wilson1D(sg[:, :, 1])
    @test isapprox(mod(2z, 1), 0.0; atol=1e-6) || isapprox(mod(2z, 1), 1.0; atol=1e-6)
end
