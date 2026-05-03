using Test
using LatticeQM
using LatticeQM.Superconductivity: BdGOperator, addpairing!
using LinearAlgebra: eigvals

# BdG construction: Nambu doubling. Starting from an N-orbital normal-state
# Hops, BdGOperator(h) lives in 2N-orbital Nambu space and the bare diagonal is
# diag(h, -conj(h)). With no pairing, eigenvalues come in ±E pairs.
@testset "Superconductivity: BdGOperator Nambu doubling and ±E spectrum" begin
    lat = Geometries.honeycomb()
    h = Hops()
    Operators.nearestneighbor!(h, lat, -1.0)
    Hbdg = BdGOperator(h)

    @test hopdim(Hbdg) == hopdim(h)               # custom: counts electron sector
    # Underlying Hops has Nambu dimension 2 × hopdim(h).
    @test size(first(values(Hbdg.h)), 1) == 2 * hopdim(h)

    # Eigenvalues at any k must be symmetric around zero (particle-hole).
    # k must match latticedim(honeycomb)=2.
    for k in ([0.1, 0.2], [0.3, -0.1], [0.0, 0.0])
        Hk = Matrix(Hbdg.h(k))
        ev = sort(real.(eigvals(Hk)))
        # Spectrum is reflection-symmetric: sort(ε) ≈ -reverse(sort(ε))
        @test maximum(abs, ev .+ reverse(ev)) < 1e-10
    end
end

# Adding a constant onsite s-wave pairing Δ opens a quasiparticle gap of size
# 2Δ at zero filling for graphene at the Dirac points (where the normal-state
# band touches zero). Test the gap size at K = (1/3, -1/3).
@testset "Superconductivity: s-wave pairing opens 2Δ gap at Dirac point" begin
    lat = Geometries.honeycomb()
    h = Hops()
    Operators.nearestneighbor!(h, lat, -1.0)
    Hbdg = BdGOperator(h)

    Δ_val = 0.25
    Δ = Hops([0, 0] => ComplexF64[Δ_val 0.0; 0.0 Δ_val])
    addpairing!(Hbdg, Δ)

    # At Dirac point the normal spectrum has E_± = 0; with constant Δ added,
    # the BdG spectrum at the Dirac point should contain ±Δ pairs.
    Hk = Matrix(Hbdg.h([1/3, -1/3]))
    evs = sort(real.(eigvals(Hk)))
    # Particle-hole symmetry preserved.
    @test maximum(abs, evs .+ reverse(evs)) < 1e-10
    # The two innermost levels are ±Δ.
    n = length(evs)
    @test isapprox(evs[n ÷ 2 + 1], +Δ_val; atol=1e-10)
    @test isapprox(evs[n ÷ 2],     -Δ_val; atol=1e-10)
end
