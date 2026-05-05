using Test
using LatticeQM

# Compute σ(ω) from the raw current-current correlator χ(ω) returned by
# `LinearResponse.opticalconductivity`. The package's low-level routine
# returns χ; the standard transformation
#
#     σ(ω) = -(χ(ω) - χ(0)) / (ω + iΓ)
#
# subtracts the diamagnetic / static piece and divides out the photon
# energy. Mirrors the recipe in `extra/examples/graphene/opticalconductivity.jl`.
σ_from_χ(χ, ωs, Γ) = -(χ .- χ[begin]) ./ (ωs .+ 1im * Γ)

# Regression: graphene NN at half-filling, μ=0, on a coarse k-grid. The Kubo
# routine + getcurrentoperators(::Hops) interplay must satisfy three exact
# (or near-exact) symmetry-protected invariants:
#
#   1. σ_xx = σ_yy. With both currents evaluated on the same k-grid this is
#      not just a C3-symmetry statement; it follows by direct identity of
#      the diagonal kernels that go into `kubo!`. Tolerance is machine-
#      precision because there's no statistical sampling difference.
#   2. Re σ_xy = 0. Time-reversal symmetric H (real graphene NN hopping) ⇒
#      vanishing transverse response at any ω. Catches sign-flips in the
#      `J_x J_y` matrix-element construction that would otherwise be silent.
#   3. Re σ_xx ≥ 0 above broadening. Causality / passivity: the dissipative
#      part of the conductivity is non-negative for ω > 0.
#
# Together these pin: getcurrentoperators(Hops) (already covered for
# Hermiticity in test_operators.jl), kubo! (matrix-element evaluation +
# Fermi-Dirac weighting), and the kspace_reduce!-driven loop in
# `opticalconductivity`.
@testset "LinearResponse: graphene σ_xx symmetries" begin
    lat = Geometries.honeycomb()
    H = Operators.graphene(lat; format=:dense, mode=:nospin)

    ωs = collect(LinRange(0.05, 3.0, 30))
    Γ, T, μ = 0.05, 0.01, 0.0

    χxx = LinearResponse.opticalconductivity(ωs, 1, 1, H, lat;
                                               klin=30, μ=μ, Γ=Γ, T=T, hidebar=true)
    χyy = LinearResponse.opticalconductivity(ωs, 2, 2, H, lat;
                                               klin=30, μ=μ, Γ=Γ, T=T, hidebar=true)
    χxy = LinearResponse.opticalconductivity(ωs, 1, 2, H, lat;
                                               klin=30, μ=μ, Γ=Γ, T=T, hidebar=true)

    σxx = σ_from_χ(χxx, ωs, Γ)
    σyy = σ_from_χ(χyy, ωs, Γ)
    σxy = σ_from_χ(χxy, ωs, Γ)

    # 1. Diagonal isotropy on the same k-grid.
    @test maximum(abs, real(σxx) .- real(σyy)) < 1e-10

    # 2. Vanishing Hall response under time reversal.
    @test maximum(abs, real(σxy)) < 1e-10

    # 3. Positive dissipative diagonal response above broadening.
    @test all(real(σxx) .> -1e-6)
end
