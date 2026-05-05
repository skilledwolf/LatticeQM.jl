using Test
using LatticeQM

# Fixture: analytic graphene conduction band ε_+(k)
include("analyticbands.jl")
using .analyticbands: grapheneconductionband

# Regression: nearest-neighbour graphene tight-binding bands must reproduce the
# analytic two-band result ε_±(k) = ±|t|·|f(k)| to numerical precision along
# the high-symmetry k-path. This pins down the joint behaviour of:
#   - Geometries.honeycomb()
#   - Operators.nearestneighbor!
#   - Operators.setfilling!
#   - kpath / Spectrum.getbands
@testset "Spectrum: graphene NN bands match analytic formula" begin
    lat = Geometries.honeycomb()
    T = Hops()
    Operators.nearestneighbor!(T, lat)
    Operators.setfilling!(T, 0.5; nk=100)

    ks = kpath(lat; num_points=200)
    bands = getbands(T, ks)
    upper_analytic = map(grapheneconductionband, eachcol(ks.points))
    lower_analytic = -upper_analytic

    mae_lower = sum(abs, bands.bands[1, :] .- lower_analytic)
    mae_upper = sum(abs, bands.bands[2, :] .- upper_analytic)

    @test mae_lower + mae_upper < 1e-12
end

# Regression: graphene NN DOS over its [-3, 3] bandwidth must satisfy three
# qualitative invariants that catch most kinds of breakage in `getdos`,
# `dos_compute!`, the `Parallel.kspace_reduce!` accumulation pipeline, and
# the eigenvalue-only `geteigvals!` shortcut they sit on top of:
#   1. Particle-hole symmetric (DOS(ω) ≈ DOS(-ω) on a bipartite lattice).
#   2. Total spectral weight = N_bands · π = 2π (the package's `dos!` accumulates
#      `Im[1/(ω − iΓ − ϵ)]` without dividing by π, so the per-band integral is
#      π rather than 1 — this is the existing normalisation convention; the
#      test asserts it explicitly so a future divide-by-π change is caught).
#   3. Van Hove peaks at ω = ±1 dominate the Dirac dip at ω = 0.
@testset "Spectrum: graphene NN DOS shape" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)

    ωs, dos = Spectrum.getdos(H, -3.05, 3.05, 401; klin=120, Γ=0.05)
    dos_real = real.(dos)

    # 1. Particle-hole symmetry to numerical precision (sublattice-symmetric H).
    @test maximum(abs, dos_real .- reverse(dos_real)) / maximum(abs, dos_real) < 1e-10

    # 2. Spectral weight (Riemann sum) ≈ 2π within ~5% — see header.
    integral = sum(dos_real) * (ωs[2] - ωs[1])
    @test isapprox(integral, 2π; rtol=0.05)

    # 3. Van Hove peaks at |ω|=1 dominate the Dirac dip at ω=0.
    i_dip = argmin(abs.(ωs .- 0.0))
    i_vH⁺ = argmin(abs.(ωs .- 1.0))
    i_vH⁻ = argmin(abs.(ωs .+ 1.0))
    @test dos_real[i_vH⁺] > 3 * dos_real[i_dip]
    @test dos_real[i_vH⁻] > 3 * dos_real[i_dip]
end
