using Test
using LatticeQM
using LinearAlgebra: dot, norm, diag

# analyticbands.jl is also included by test_spectrum_graphene.jl; redefining
# the module is harmless but using-imports must avoid Main-level conflicts.
include("analyticbands.jl")
import .analyticbands: grapheneconductionband

# Hartree-Fock self-consistency for graphene Hubbard at U/t = 4.
# This is a slow integration test (~30-60 s); set LATTICEQM_SKIP_SLOW=1 to skip.
function get_graphene_HMF(U=4.0; filling=0.5, init=:antiferro, T=0.01, β=0.20)
    lat = Geometries.honeycomb()
    hops = Operators.graphene(lat; mode=:spinhalf)

    v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=U)

    ρ_init = Meanfield.initialguess(v, init; lat=lat)

    # NOTE: clear_trace=false (used in legacy hubbardtest.jl) is not accepted
    # by the current fixedpoint! signature.
    ρ, _, HMF, converged, _ = Meanfield.solvehartreefock(
        hops, v, ρ_init, filling; klin=30, iterations=800, tol=1e-5,
        T=T, β=β, show_trace=false)

    return lat, ρ, HMF, converged
end

# Sanity: at vanishing U, the HF self-consistent bands must reproduce the
# non-interacting graphene bands. We use U = 1e-6 (not exactly zero) because
# Operators.gethubbard internally drops entries below sqrt(eps()) ≈ 1.5e-8;
# U = 0 would yield an empty Hops and break Meanfield.initialguess.
function meanfieldsanitycheck()
    lat, _, HMF, converged = get_graphene_HMF(1e-6; T=0.0, β=0.70)
    converged || return false

    ks = kpath(lat; num_points=200)
    nk = size(ks)[2]
    bands_mf = getbands(HMF.h, ks)
    bands_mf.bands .-= HMF.μ

    upper = map(grapheneconductionband, eachcol(ks.points))
    lower = -upper

    for band in eachrow(bands_mf.bands)
        mae_upper = sum(abs, band .- upper) / nk
        mae_lower = sum(abs, band .- lower) / nk
        # Modest tolerance: density-matrix arithmetic introduces fp noise.
        if !(mae_upper < 6e-4) && !(mae_lower < 6e-4)
            return false
        end
    end
    return true
end

function get_gap_at_U(U=4.0; filling=0.5, init=:antiferro, T=0.01, β=0.20)
    lat, _, HMF, _ = get_graphene_HMF(U; filling=filling, init=init, T=T, β=β)
    ks = kpath(lat; num_points=200)
    return Spectrum.bandgap_filling(HMF.h, ks, filling)
end

# The Hubbard model on graphene has a Mott transition near U_c ≈ 2.2 t.
# The location of the maximum curvature of the gap-vs-U curve is a robust
# numerical fingerprint of the SCF solver. Start above 0 because
# Operators.gethubbard returns an empty Hops at U=0 (entries below
# sqrt(eps()) are dropped), which breaks Meanfield.initialguess.
function hubbardmodelcriticalpoint()
    Us = range(1e-6; stop=4.0, length=10)
    _, minidx = findmin(U -> abs(U - 2.23), Us)
    gaps = map(U -> get_gap_at_U(U; init=:antiferro), Us)

    dgaps = diff(gaps)
    ddgaps = diff(dgaps)
    maxdd = maximum(ddgaps)

    return ddgaps[minidx - 2] == maxdd || ddgaps[minidx - 1] == maxdd
end

# Anderson acceleration must converge to the same ground-state energy as
# linear mixing on graphene Hubbard at U=4, but in far fewer iterations. We
# don't pin the exact iteration count (it depends on system, history depth,
# damping), only that Anderson stays at least 3× faster than the β=0.20
# linear baseline that the Néel test uses.
function andersonacceleration()
    lat = Geometries.honeycomb()
    hops = Operators.graphene(lat; mode=:spinhalf)
    v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=4.0)
    ρ_init = Meanfield.initialguess(v, :antiferro; lat=lat)

    log_lin = Tuple{Int,Float64,Float64}[]
    _, ϵ_lin, _, conv_lin, _ = Meanfield.solvehartreefock(
        hops, v, deepcopy(ρ_init), 0.5; klin=15, iterations=400, tol=1e-7,
        T=0.01, β=0.20, verbose=false, show_trace=false,
        log_callback=(it, ϵ, r) -> push!(log_lin, (it, ϵ, r)))

    log_and = Tuple{Int,Float64,Float64}[]
    _, ϵ_and, _, conv_and, _ = Meanfield.solvehartreefock(
        hops, v, deepcopy(ρ_init), 0.5; klin=15, iterations=400, tol=1e-7,
        T=0.01, β=1.0, acceleration=:anderson, anderson_depth=5, verbose=false,
        show_trace=false,
        log_callback=(it, ϵ, r) -> push!(log_and, (it, ϵ, r)))

    converged = conv_lin && conv_and
    same_energy = isapprox(ϵ_lin, ϵ_and; atol=1e-6)
    much_faster = length(log_and) * 3 ≤ length(log_lin)
    return converged && same_energy && much_faster
end

# Regression test for the purification SCF. Until commit 6ffeb8f-ish, the
# closure inside `_solveselfconsistent_purification_impl!` did
# `ρ1 = canonicalpurification_grid_realspace(...)`, rebinding the local
# parameter. fixedpoint!'s outer ρ1 was never updated, the residual was 0
# from iter 1, and the routine returned the initial guess as the "converged"
# answer.
#
# We cover that bug here without committing to numerical agreement with the
# diagonalization solver — the canonical-purification algorithm in the tree
# has further unrelated issues (notably `compute_AB_product` truncates the
# real-space convolution to keys(A), and the inner McWeeny `c`-coefficient
# can fall outside [0,1] on some Hamiltonians) that make a hard energy-match
# test infeasible until that algorithm is fixed.
function purificationiterates()
    lat = Geometries.honeycomb_twisted(5)
    hops = Operators.graphene(lat; format=:sparse, mode=:spinhalf, tz=0.52)
    v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=2.0)
    ρ_init = Meanfield.initialguess(v, :ferro; lat=lat)

    log = Tuple{Int,Float64,Float64}[]
    ρ, ϵ, _, _, _ = Meanfield.solvehartreefock_purification(
        hops, v, deepcopy(ρ_init), 0.5; klin=3, iterations=6, tol=1e-2,
        β=0.5, relative=false, verbose=false,
        log_callback=(it, ϵ, r) -> push!(log, (it, ϵ, r)))

    energies = [step[2] for step in log]
    residuals = [step[3] for step in log]

    iterated = length(log) ≥ 2 &&
               length(unique(round.(energies, sigdigits=6))) ≥ 2 &&
               length(unique(round.(residuals, sigdigits=6))) ≥ 2
    zk = LatticeQM.TightBinding.zerokey(ρ)
    finite = isfinite(real(ϵ)) && all(isfinite, real.(diag(ρ[zk])))
    energy_decreased = real(energies[end]) < real(energies[1])
    hermitian = LatticeQM.TightBinding.ishermitian(ρ)

    return iterated && finite && energy_decreased && hermitian
end

# At U=4, the AF Hartree-Fock solution should be a Néel state: the magnetisation
# vectors on the two sublattices must be exactly anti-aligned.
function hubbardmodelantiferro()
    lat, ρ, _, _ = get_graphene_HMF()
    suba = getoperator(lat, "sublatticeAspin")
    subb = getoperator(lat, "sublatticeBspin")

    subamag_ops = [suba * s for s in getoperator(lat, ["sx", "sy", "sz"])]
    subbmag_ops = [subb * s for s in getoperator(lat, ["sx", "sy", "sz"])]
    subamag = real([Operators.trace(ρ, mag) for mag in subamag_ops])
    subbmag = real([Operators.trace(ρ, mag) for mag in subbmag_ops])

    # Tolerance 1e-8 matches the residual SCF asymmetry left at tol=1e-5
    # (legacy code asked for 1e-12 which is tighter than the SCF residual).
    return isapprox(dot(subamag, subbmag) / norm(subamag) / norm(subbmag),
                    -1.0; atol=1e-8)
end

@testset "Meanfield: Hubbard graphene (slow, opt-in)" begin
    # Default Pkg.test() skips this set: it's a 5-minute SCF integration test.
    # Opt in with LATTICEQM_FULL_TESTS=1 (CI nightly job, release candidate
    # checks). The fast suite covers the Meanfield API surface indirectly via
    # the analytic-band tests; this one is the end-to-end SCF regression.
    if get(ENV, "LATTICEQM_FULL_TESTS", "0") != "1"
        @info "Skipping slow Hubbard mean-field tests (set LATTICEQM_FULL_TESTS=1 to run)"
    else
        @testset "U=0 reproduces non-interacting bands" begin
            @test meanfieldsanitycheck()
        end
        @testset "Mott critical-point fingerprint near U≈2.2" begin
            @test hubbardmodelcriticalpoint()
        end
        @testset "Néel order at U=4 (sublattice mags anti-aligned)" begin
            @test hubbardmodelantiferro()
        end
        @testset "Anderson acceleration converges 3× faster than β=0.20 linear" begin
            @test andersonacceleration()
        end
        @testset "Purification SCF actually iterates (regression for in-place fix)" begin
            @test purificationiterates()
        end
    end
end
