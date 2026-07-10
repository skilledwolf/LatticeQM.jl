using Test
using LatticeQM
using LinearAlgebra: I, tr, Diagonal, Hermitian, eigen, eigvals, dot

# Spin operators on a 2-orbital honeycomb sit in a (2N)×(2N) = 4×4 space.
# They must satisfy the spin-1/2 Pauli algebra block-by-block.
@testset "Operators: spin operators on honeycomb" begin
    lat = Geometries.honeycomb()
    N = 2
    sx = getoperator(lat, "sx")
    sy = getoperator(lat, "sy")
    sz = getoperator(lat, "sz")

    @testset "shape" begin
        @test size(sx) == (2N, 2N)
        @test size(sy) == (2N, 2N)
        @test size(sz) == (2N, 2N)
    end

    @testset "Hermiticity" begin
        @test sx ≈ sx'
        @test sy ≈ sy'
        @test sz ≈ sz'
    end

    @testset "involution: sᵢ² = I (in the 2N×2N space)" begin
        Imat = Matrix{ComplexF64}(I, 2N, 2N)
        @test sx * sx ≈ Imat
        @test sy * sy ≈ Imat
        @test sz * sz ≈ Imat
    end

    @testset "commutation: [sx, sy] = 2i sz (and cyclic)" begin
        @test sx * sy - sy * sx ≈ 2im * sz
        @test sy * sz - sz * sy ≈ 2im * sx
        @test sz * sx - sx * sz ≈ 2im * sy
    end

    @testset "name-aliasing: 'SX' / 'sx' / 'MX' all return same operator" begin
        @test getoperator(lat, "SX") ≈ sx
        @test getoperator(lat, "MX") ≈ sx
    end

    @testset "spin up/down projectors sum to identity and are orthogonal" begin
        sup = getoperator(lat, "spinup")
        sdn = getoperator(lat, "spindown")
        Imat = Matrix{ComplexF64}(I, 2N, 2N)
        @test sup + sdn ≈ Imat
        @test sup * sdn ≈ zeros(2N, 2N)
        @test sup * sup ≈ sup       # idempotent
    end
end

# Sublattice projectors must be orthogonal idempotents that partition the
# Hilbert space.
@testset "Operators: sublattice projectors on honeycomb" begin
    lat = Geometries.honeycomb()
    PA = getoperator(lat, "sublatticeA")
    PB = getoperator(lat, "sublatticeB")
    @test size(PA) == (2, 2)        # honeycomb has N=2 orbitals; spinless here
    @test PA + PB ≈ Matrix{Float64}(I, 2, 2)
    @test PA * PB ≈ zeros(2, 2)
    @test PA * PA ≈ PA
    @test tr(PA) ≈ 1
    @test tr(PB) ≈ 1
end

# Sublattice-spin projectors live in (2N)×(2N) and project onto e.g. all
# A-sublattice spin states.
@testset "Operators: sublattice-spin projectors" begin
    lat = Geometries.honeycomb()
    PAs = getoperator(lat, "sublatticeAspin")
    PBs = getoperator(lat, "sublatticeBspin")
    @test size(PAs) == (4, 4)
    @test PAs + PBs ≈ Matrix{Float64}(I, 4, 4)
    @test PAs * PBs ≈ zeros(4, 4)
    @test tr(PAs) ≈ 2     # 1 sublattice site × 2 spin states
end

# trace and expval on Hops; uses the onsite (zero-key) sector.
@testset "Operators: trace and expval on identity Hops" begin
    Imat = Matrix{ComplexF64}(I, 4, 4)
    H = Hops([0, 0] => Imat)
    @test Operators.trace(H) ≈ 4
    @test Operators.expval(H) ≈ tr(Imat)   # sum of all entries — onsite only here
end

# Regression: expval must compute Σ_L tr[ρ_L O_{−L}] = (1/Nk) Σ_k tr[ρ(k)O(k)].
# The old L↔L pairing computed tr[ρ_L O_Lᵀ], which (a) flips the sign of every
# transpose-odd onsite operator (sy!) and (b) is wrong for all inter-cell blocks.
@testset "Operators: expval matches direct Bloch contraction (incl. sy)" begin
    lat = Geometries.honeycomb()
    H = Operators.graphene(lat; mode=:spinhalf)
    Operators.addzeeman!(H, lat, [0.0, 0.4, 0.0])

    nk = 9
    ks = hcat([[i / nk, j / nk] for i in 0:nk-1, j in 0:nk-1][:]...)
    μ = -0.5   # inside the lower bands: nonzero polarization, no state at ϵ=μ

    # getdensitymatrix! accumulates only into pre-existing keys.
    ρ = Hops(Dict(L => zeros(ComplexF64, size(H[L])) for L in keys(H)))
    Operators.getdensitymatrix!(ρ, H, ks, μ; T=0.0)

    # Independent reference: plain Bloch-sum expectation values.
    Sy = getoperator(lat, "sy")
    sy_direct = 0.0
    eband_direct = 0.0
    for k in eachcol(ks)
        hk = Matrix(H(k))
        ϵs, U = eigen(Hermitian(hk))
        ρk = U * Diagonal(float.(ϵs .< μ)) * U'
        sy_direct += real(tr(ρk * Sy))
        eband_direct += real(tr(ρk * hk))
    end
    sy_direct /= size(ks, 2)
    eband_direct /= size(ks, 2)

    @test abs(sy_direct) > 1e-3                                   # test has teeth
    @test real(Operators.expval(ρ, "sy", lat)) ≈ sy_direct atol=1e-10
    # Inter-cell pairing: contracting ρ against the Hamiltonian itself gives
    # the occupied band energy only if L pairs with −L.
    @test real(Operators.expval(ρ, H)) ≈ eband_direct atol=1e-10
    # magnetization wraps expval for all three components.
    mx, my, mz = real.(Operators.magnetization(ρ, lat))
    @test my ≈ sy_direct atol=1e-10
    @test abs(mx) < 1e-8
    @test abs(mz) < 1e-8
end

# Regression: Spectrum/Hops current operators J_α(k) = ∂H(k)/∂k_α must be
# Hermitian for any Hermitian H. The Hops-input branch of getcurrentoperators
# is the documented "preferred path" for tight-binding Hamiltonians, so this
# is the right one to pin. Catches sign-flip / position-convention drift in
# `currentoperators!` and the dense/sparse helpers it dispatches to.
@testset "Operators: getcurrentoperators(Hops) is Hermitian" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    Js = Operators.getcurrentoperators(lat, H)

    # One operator per lattice direction (graphene = 2). The function
    # previously returned `spacedim(lat) = 3` operators with the third one
    # identically zero; current.jl now returns only the meaningful ones.
    @test length(Js) == LatticeQM.Structure.Lattices.latticedim(lat)

    for k in ([0.0, 0.0], [1/3, -1/3], [0.5, 0.5], [-0.2, 0.7])
        for J in Js
            Jk = Matrix(J(k))
            @test maximum(abs, Jk .- Jk') < 1e-10
        end
    end
end

# ---------------------------------------------------------------------------
# Peierls substitution.
# ---------------------------------------------------------------------------

# All Peierls entry points share one sign convention:
#   t_ij → t_ij · exp(+i 2π/φ0 ∫_{r_j}^{r_i} A·dl).
# Pin both gauges against a brute-force numerical line integral.
@testset "Operators: uniformfieldphase = +∫A·dl (both gauges)" begin
    # symmetric gauge A = ½ B×r
    B3 = [0.3, -0.7, 0.9]
    r1 = [0.8, -0.4, 0.25]
    r2 = [-0.3, 0.9, -0.6]
    A_sym(r) = 0.5 .* [B3[2]*r[3]-B3[3]*r[2], B3[3]*r[1]-B3[1]*r[3], B3[1]*r[2]-B3[2]*r[1]]
    nseg = 20_000
    δ = r1 .- r2
    lineint = sum(dot(A_sym(r2 .+ (t - 0.5) / nseg .* δ), δ ./ nseg) for t in 1:nseg)
    @test Operators.uniformfieldphase(r1, r2; B=B3) ≈ lineint atol=1e-8

    # 2D positions: treated as z=0 plane (this path used to throw a
    # DimensionMismatch and, for 3D input, returned −∫A·dl).
    @test Operators.uniformfieldphase([1.0, 0.0], [0.0, 1.0]; B=[0, 0, 1.0]) ≈
          Operators.uniformfieldphase([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]; B=[0, 0, 1.0])

    # in-plane gauge A = (B2 z, −B1 z, 0)
    B2 = [0.4, -1.1]
    A_ip(r) = [B2[2]*r[3], -B2[1]*r[3], 0.0]
    lineint_ip = sum(dot(A_ip(r2 .+ (t - 0.5) / nseg .* δ), δ ./ nseg) for t in 1:nseg)
    @test Operators.uniformfieldphase_inplane(r1, r2; B=B2) ≈ lineint_ip atol=1e-8
end

# A uniform field in a fixed gauge breaks lattice translation invariance;
# the old code silently produced a non-Hermitian Bloch Hamiltonian (which the
# dense solver then silently symmetrized into wrong spectra). Now it refuses.
@testset "Operators: peierls! rejects uniform B on periodic lattices" begin
    lat = Geometries.honeycomb()
    H = Hops()
    Operators.nearestneighbor!(H, lat, -1.0)
    @test_throws ArgumentError Operators.peierls!(H, lat, [0.0, 0.0, 0.3])
end

# On a finite (0D) lattice the symmetric gauge is valid. Check Hermiticity,
# the gauge-invariant plaquette flux (pins the +∫A·dl sign), and the exact
# 4-site-ring spectrum.
@testset "Operators: peierls! on a 0D flake — flux and spectrum" begin
    Φ = 0.3
    # unit square, counterclockwise: 1=(0,0), 2=(1,0), 3=(1,1), 4=(0,1)
    pos = Float64[0 1 1 0; 0 0 1 1]
    lat = LatticeQM.Structure.Lattices.Lattice(Matrix(1.0I, 2, 2), 0, pos, zeros(0, 4))

    ring = zeros(ComplexF64, 4, 4)
    for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 1))
        ring[i, j] = -1.0
        ring[j, i] = -1.0
    end
    H = Hops(Int[] => copy(ring))
    Operators.peierls!(H, lat, [0.0, 0.0, Φ])
    M = H[Int[]]

    @test M ≈ M'                                   # Hermitian (was broken globally)

    # counterclockwise plaquette product: t⁴·exp(+i2π Φ_enclosed), t = −1
    P = M[2, 1] * M[3, 2] * M[4, 3] * M[1, 4]
    @test angle(P) ≈ 2π * Φ atol=1e-12

    # gauge-invariant: spectrum equals the uniform-phase 4-ring result
    ϵ = eigvals(Hermitian(M))
    ϵ_exact = sort([-2 * cos(2π * (n + Φ) / 4) for n in 0:3])
    @test ϵ ≈ ϵ_exact atol=1e-12
end

# ---------------------------------------------------------------------------
# Rhombohedral (ABC) trilayer Slonczewski-Weiss-McClure hoppings.
# Regression: γ3 used to require in-plane distance 2a, which matches no
# physical A1-B2 bond (they sit at in-plane distance a, like γ4) — the term
# was silently absent and a spurious far shell got γ3 instead.
# ---------------------------------------------------------------------------
@testset "Operators: graphene_rhombohedral puts γ3/γ4 on the right bonds" begin
    a, d = 1.0, 3.0
    γ0, γ1, γ2, γ3, γ4 = -3.16, 0.502, -0.0171, -0.377, -0.099
    lat = Geometries.honeycomb_ABC(a, d)
    H = Operators.graphene_rhombohedral(lat; a=a, d=d, γ0=γ0, γ1=γ1, γ2=γ2, γ3=γ3, γ4=γ4)

    A = LatticeQM.Structure.Lattices.getA(lat)
    R = LatticeQM.Structure.Lattices.allpositions(lat)
    N = LatticeQM.Structure.Lattices.countorbitals(lat)

    # classify every stored bond by geometry and check its amplitude
    seen_γ3 = 0
    seen_γ4 = 0
    for L in keys(H), i in 1:N, j in 1:N
        amp = H[L][i, j]
        abs(amp) < 1e-12 && continue
        δr = (R[1:3, i] .+ A * L) .- R[1:3, j]
        rin = LinearAlgebra.norm(δr[1:2])
        dz = abs(δr[3])
        δsub = abs(R[4, i] - R[4, j])
        if isapprox(dz, d; atol=1e-8) && isapprox(rin, a; atol=1e-8)
            if δsub > 0.5
                @test real(amp) ≈ γ3 atol=1e-12   # A1-B2 family
                seen_γ3 += 1
            else
                @test real(amp) ≈ γ4 atol=1e-12   # A1-A2 / B1-B2 family
                seen_γ4 += 1
            end
        elseif isapprox(rin, 2a; atol=1e-8)
            @test false   # nothing lives on the old spurious 2a shell
        end
    end
    # each of the 4 non-dimer interlayer interfaces contributes 3 bonds ×
    # both hopping directions; just require the terms exist in numbers
    @test seen_γ3 >= 12
    @test seen_γ4 >= 12

    @test LatticeQM.TightBinding.ishermitian(H)
end

# Haldane term must not silently vanish when the bond length is not 1.
@testset "Operators: Haldane term is bond-length aware" begin
    t2 = 0.1
    H1 = Operators.gethaldane(Geometries.honeycomb(), t2)
    H2 = Operators.gethaldane(Geometries.honeycomb(1.42), t2)   # a ≠ 1
    n1 = sum(count(!iszero, M) for M in values(H1.data))
    n2 = sum(count(!iszero, M) for M in values(H2.data))
    @test n1 > 0
    @test n2 == n1        # same connectivity, rescaled geometry
end
