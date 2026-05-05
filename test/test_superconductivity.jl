using Test
using LatticeQM
using LatticeQM.Superconductivity: BdGOperator, addpairing!
using LinearAlgebra: eigvals, tr, I, Hermitian
import LatticeQM.TightBinding

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

# Validation regressions for the construction-time guards in BdGOperator /
# addpairing!. Each of these should error with a clear message rather than
# silently constructing a non-Hermitian Nambu Hamiltonian.
@testset "Superconductivity: validation guards" begin
    lat = Geometries.honeycomb()
    h = Hops(); Operators.nearestneighbor!(h, lat, -1.0)

    @testset "non-Hermitian h is rejected" begin
        h_bad = Hops()
        # Asymmetric: only the +1 offset, no −1 partner.
        h_bad[[1, 0]]  = ComplexF64[0.0 1.0; 0.0 0.0]
        h_bad[[0, 0]]  = ComplexF64[0.0 0.0; 0.0 0.0]
        @test_throws AssertionError BdGOperator(h_bad)
    end

    @testset "pairing dict missing -R partner is rejected" begin
        H = BdGOperator(h)
        Δ_oneside = Hops()
        Δ_oneside[[1, 0]] = ComplexF64[0.1 0.0; 0.0 0.1]
        # Note: no Δ_oneside[[-1, 0]] entry → addpairing! must error.
        @test_throws ErrorException addpairing!(H, Δ_oneside)
    end

    @testset "Δ[+R] and Δ[-R] independent: BdG still Hermitian" begin
        # The construction `[0 Δ[R]; Δ[-R]' 0]` produces a Hermitian BdG
        # *regardless* of any Δ[+R] vs Δ[-R] relationship — they parametrise
        # independent upper-right blocks. addpairing! must NOT reject this
        # (the SCF's natural ΔMF lives here).
        H = BdGOperator(h)
        Δ_indep = Hops()
        Δ_indep[[1, 0]]  = ComplexF64[0.1 0.2; 0.0 0.1]
        Δ_indep[[-1, 0]] = ComplexF64[0.1 0.0; 0.0 0.1]   # NOT Δ[+R]'
        addpairing!(H, Δ_indep)
        @test TightBinding.ishermitian(H.h)
    end
end

# Bond pairing exercise: a NN-bond Hermitian-compatible Δ on graphene.
# Rather than the trivial onsite-identity test, this actually walks the
# offset machinery in addpairing! (R and −R both present) and verifies the
# resulting BdG is Hermitian and particle-hole symmetric.
@testset "Superconductivity: bond pairing — Hermiticity and ±E spectrum" begin
    lat = Geometries.honeycomb()
    h = Hops(); Operators.nearestneighbor!(h, lat, -1.0)

    # Constant bond pairing on (+1, 0) and its partner (−1, 0). The values
    # satisfy Δ[−R] = Δ[R]† by construction.
    Δ_val = 0.05 + 0.0im
    bond = ComplexF64[0.0 Δ_val; 0.0 0.0]
    Δ = Hops()
    Δ[[ 1,  0]] = bond
    Δ[[-1,  0]] = bond'   # Δ[-R]' = bond — required by addpairing!

    Hbdg = BdGOperator(h, Δ)
    @test TightBinding.ishermitian(Hbdg.h)

    for k in ([0.1, 0.2], [0.3, -0.1], [1/3, -1/3])
        Hk = Matrix(Hbdg.h(k))
        @test maximum(abs, Hk .- Hk') < 1e-12   # Hermitian at every k
        ev = sort(real.(eigvals(Hk)))
        @test maximum(abs, ev .+ reverse(ev)) < 1e-10  # ±E pairs
    end
end

# Energy convention: the BdG `getdensitymatrix!` returns the *physical*
# grand-canonical kinetic ⟨H_phys − μN̂⟩, recovered from the BdG sum via the
# Nambu c-number correction. Compare to a direct many-body computation on a
# tiny model where both sides can be evaluated exactly.
@testset "Superconductivity: BdG kinetic energy matches direct ⟨H_phys⟩" begin
    # Two-orbital, single k-point (Γ), no pairing — easiest case to compute
    # by hand and compare to the BdG-corrected output.
    lat = Geometries.honeycomb()
    h = Hops(); Operators.nearestneighbor!(h, lat, -1.0)
    Hbdg = BdGOperator(h)

    # Single k-point at Γ. h(0) = nearest-neighbour graphene at Γ — bands at
    # ±3 with eigenvectors (1, ±1)/√2 in the AB sublattice basis.
    ks = reshape([0.0, 0.0], 2, 1)

    # Compute via getdensitymatrix! at μ=0, T→0 (use small T for numerics).
    ρ = BdGOperator(zero(LatticeQM.Utils.copyelectronsector(Hbdg)))
    ϵ_kin = Operators.getdensitymatrix!(ρ, Hbdg, ks, 0.0;
                                        T=1e-6, hidebar=true, format=:dense,
                                        multimode=:serial)

    # Direct: at μ=0 and T→0, the lower band (E = −3) is filled, upper is
    # empty. ⟨H⟩ = E_filled = −3.
    @test isapprox(ϵ_kin, -3.0; atol=1e-3)
end

# Filling regression: verify Spectrum.filling(::BdGOperator) tracks the
# physical electron count even when Δ ≠ 0 — this is the path the SCF takes
# to fit μ self-consistently with the gap turned on. With no pairing the
# answer is the normal-state filling at μ; with weak pairing the deviation
# is O(Δ²/W), still close to the target.
@testset "Superconductivity: BdG filling tracks target with weak Δ" begin
    lat = Geometries.honeycomb()
    h = Hops(); Operators.nearestneighbor!(h, lat, -1.0)
    Hbdg = BdGOperator(h)

    Δ_val = 0.05
    Δ = Hops([0, 0] => ComplexF64[Δ_val 0.0; 0.0 Δ_val])
    addpairing!(Hbdg, Δ)

    # Use a small but representative grid; at T = 0.05 the Fermi function
    # is smooth enough that a 12² grid resolves the half-filling point.
    ks = LatticeQM.Structure.regulargrid(nk=12^2)
    f0 = Spectrum.filling(Hbdg, LatticeQM.Structure.points(ks), 0.0; T=0.05)
    # At μ=0 with particle-hole symmetric h, filling should be exactly 1/2
    # to within Δ²/W corrections (Δ=0.05, W≈3 → ~3×10⁻⁴).
    @test isapprox(f0, 0.5; atol=5e-3)
end
