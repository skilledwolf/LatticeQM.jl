using Test
using LatticeQM
using LatticeQM.Floquet:
    periodicDrive, addFreq!, addcos!, getFloquetMatrix, getfirstFBZ, transform, transformeigvecs
using LinearAlgebra: eigvals, eigen, Hermitian, I
import LinearAlgebra  # qualify ishermitian: the bare name clashes with LatticeQM's Hops method when the whole suite shares Main

# Smoke + sanity for Floquet construction. With no harmonics, the Floquet
# Hamiltonian for a static H₀ over (2M+1) Fourier replicas is block-diagonal
# with H₀ + n·ω·I on the n-th block, so its eigenvalues are
#   { εᵢ + n·ω : i = 1..N, n = -M..M }.
@testset "Floquet: empty drive ⇒ block-diagonal n·ω-shifted spectrum" begin
    ω = 1.7
    H0 = ComplexF64[1.0 0.5+0.2im; 0.5-0.2im -0.4]
    drive = periodicDrive(ω)              # no harmonics added
    @test drive.omega == ω
    @test isempty(drive.ns)
    @test isempty(drive.operators)

    M = 2
    HF = getFloquetMatrix(M, H0, drive)
    static_evs = real.(eigvals(H0))
    expected = sort([ε + n * ω for ε in static_evs for n in -M:M])
    got = sort(real.(eigvals(Matrix(HF))))
    @test maximum(abs, got .- expected) < 1e-10
end

@testset "Floquet: addFreq! and addcos! mutate harmonic content" begin
    drive = periodicDrive(2.0)
    op = ComplexF64[0 1; 1 0]
    addFreq!(drive, op, 1)
    @test drive.ns == [1]
    @test length(drive.operators) == 1

    addcos!(drive, op, 2)                # appends two harmonics: ±2
    @test sort(drive.ns) == [-2, 1, 2]
    @test length(drive.operators) == 3
end

# V(t) = V·cos(ωt) has Fourier components V_{±1} = V/2. The Floquet matrix must
# be exactly Hermitian (the dense eigensolver path wraps it in Hermitian(...),
# which reads only the upper triangle), and storing only one member of a ±n
# pair must be completed automatically via V_{-n} = V_n†.
@testset "Floquet: Hermitian completion of drive harmonics" begin
    ω = 1.9
    M = 2
    H0 = ComplexF64[0.2 0.0; 0.0 -0.2]
    V = ComplexF64[0.0 0.7-0.1im; 0.7+0.1im 0.0]   # Hermitian drive operator

    both = periodicDrive(ω)
    addcos!(both, V, 1)                            # stores both ±1 explicitly
    plusonly = periodicDrive(ω)
    addFreq!(plusonly, V / 2, 1)                   # stores only +1
    minusonly = periodicDrive(ω)
    addFreq!(minusonly, V' / 2, -1)                # stores only -1

    HFboth = Matrix(getFloquetMatrix(M, H0, both))
    HFplus = Matrix(getFloquetMatrix(M, H0, plusonly))
    HFminus = Matrix(getFloquetMatrix(M, H0, minusonly))

    @test LinearAlgebra.ishermitian(HFboth)
    @test LinearAlgebra.ishermitian(HFplus)
    @test LinearAlgebra.ishermitian(HFminus)
    @test HFplus == HFboth
    @test HFminus == HFboth

    # Regression for the Hermitian(...) upper-triangle path: a drive stored
    # only as n = -1 used to be silently discarded there.
    @test sort(eigvals(Hermitian(HFminus))) ≈ sort(eigvals(Hermitian(HFboth)))
    @test !(sort(eigvals(Hermitian(HFminus))) ≈
            sort(eigvals(Hermitian(Matrix(getFloquetMatrix(M, H0, periodicDrive(ω)))))))

    # inconsistent ±1 harmonics (non-Hermitian V(t)) are rejected
    bad = periodicDrive(ω)
    addFreq!(bad, V / 2, 1)
    addFreq!(bad, V / 2 + ComplexF64[0.3 0; 0 0], -1)
    @test_throws ArgumentError getFloquetMatrix(M, H0, bad)
end

@testset "Floquet: periodicDrive convenience constructor dispatches" begin
    ω = 2.0
    V = ComplexF64[0.0 1.0; 1.0 0.0]
    # Vector{Float64} + Matrix{ComplexF64} used to miss the convenience method
    # (invariant type parameters) and die in the default constructor.
    drive = periodicDrive(ω, [0.5, 0.0, 0.5], V)   # cos(ωt) drive, modes -1,0,1
    @test drive isa periodicDrive
    @test drive.omega == ω
    @test collect(drive.ns) == [-1, 0, 1]
    @test drive.operators[1] ≈ 0.5 * V
    @test drive.operators[3] ≈ 0.5 * V

    # ... and it produces the same Floquet matrix as adding ±1 explicitly
    M = 2
    H0 = ComplexF64[0.1 0.0; 0.0 -0.1]
    ref = periodicDrive(ω)
    addcos!(ref, V, 1)
    @test Matrix(getFloquetMatrix(M, H0, drive)) ≈ Matrix(getFloquetMatrix(M, H0, ref))
end

# Independent check: build the dense Floquet matrix for V(t) = V·cos(ωt) by
# hand — blocks H0 + pω·I on the diagonal (p = M..-M top-left to bottom-right)
# and V/2 on the first block off-diagonals — and compare quasienergies.
@testset "Floquet: quasienergies match hand-built dense Floquet matrix" begin
    ω = 2.2
    M = 3
    d = 2
    H0 = ComplexF64[0.3 0.1im; -0.1im -0.2]
    V = ComplexF64[0.0 0.4; 0.4 0.0]

    drive = periodicDrive(ω)
    addcos!(drive, V, 1)
    HF = Matrix(getFloquetMatrix(M, H0, drive))

    N = (2M + 1) * d
    Hdense = zeros(ComplexF64, N, N)
    for (bi, p) in enumerate(M:-1:-M)
        r = (bi-1)*d+1 : bi*d
        Hdense[r, r] = H0 + p * ω * I
    end
    for bi in 1:2M
        r = (bi-1)*d+1 : bi*d
        r2 = bi*d+1 : (bi+1)*d
        Hdense[r, r2] += V / 2      # couples Fourier index p to p-1
        Hdense[r2, r] += V' / 2
    end

    @test HF ≈ Hdense
    @test sort(eigvals(Hermitian(HF))) ≈ sort(eigvals(Hermitian(Hdense)))
end

@testset "Floquet: getfirstFBZ selects quasienergies in (-ω/2, ω/2]" begin
    ω = 1.0
    M = 2
    H0 = ComplexF64[0.1 0.0; 0.0 0.6]   # 0.6 crosses the FBZ edge: folds to -0.4
    HF = Matrix(getFloquetMatrix(M, H0, periodicDrive(ω)))
    evs = sort(eigvals(Hermitian(HF)))
    bands = repeat(evs, 1, 3)           # pretend three k points

    F = getfirstFBZ(bands, M, ω)
    @test size(F) == (2, 3)
    @test all(-ω/2 .< F .<= ω/2)
    @test F[:, 1] ≈ [-0.4, 0.1]

    # legacy two-argument method keeps the old middle-rows behavior
    # (documented as only valid when no band crosses the FBZ edge)
    Flegacy = getfirstFBZ(bands, M)
    @test Flegacy[:, 1] ≈ [0.1, 0.6]

    # fallback: no eigenvalue inside the window (severe truncation) — the
    # closest state is picked and folded back into the FBZ
    bands_edge = reshape([0.6, 1.6, 2.6], 3, 1)   # d = 1, M = 1
    Fe = getfirstFBZ(bands_edge, 1, ω)
    @test Fe ≈ reshape([-0.4], 1, 1)
end

# Time reconstruction: block m of a Floquet eigenvector carries the phase
# exp(-i(m-M)ωt). With ω ≠ 1 this only works if transform receives ω. For an
# undriven system every Floquet eigenpair (λ, u) must reproduce the exact
# time evolution: exp(-iλt)·transform(u, t) == exp(-i·t·H0)·ψ(0).
@testset "Floquet: transform uses ω in the reconstruction phases (ω ≠ 1)" begin
    ω = 2.5
    M = 2
    d = 2
    H0 = ComplexF64[1.0 0.5; 0.5 -0.4]
    HF = Matrix(getFloquetMatrix(M, H0, periodicDrive(ω)))
    sol = eigen(Hermitian(HF))

    for j in (1, 3, 2M*d + 2)           # eigenvectors living in different Fourier blocks
        u = sol.vectors[:, j]
        λ = sol.values[j]
        ψ0 = transform(u, d, M, ω, 0.0)
        for t in (0.3, 1.234, 2π/ω)
            ψt = exp(-im * λ * t) .* transform(u, d, M, ω, t)
            @test ψt ≈ exp(-im * t * H0) * ψ0 atol = 1e-8
        end
    end

    # driven case: the periodic part φ(t) = e^{iεt}ψ(t) must be T-periodic,
    # with T = 2π/ω — this fails if the phases miss the factor ω
    drive = periodicDrive(ω)
    addcos!(drive, ComplexF64[0.0 0.3; 0.3 0.0], 1)
    U = eigen(Hermitian(Matrix(getFloquetMatrix(M, H0, drive)))).vectors
    T = 2π / ω
    @test transformeigvecs(U, M, ω, T) ≈ transformeigvecs(U, M, ω, 0.0) atol = 1e-8
end
