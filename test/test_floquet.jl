using Test
using LatticeQM
using LatticeQM.Floquet: periodicDrive, addFreq!, addcos!, getFloquetMatrix
using LinearAlgebra: eigvals

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
