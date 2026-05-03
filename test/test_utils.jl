using Test
using LatticeQM
using LatticeQM.Utils: σ0, σ1, σ2, σ3, σX, σY, σZ, σUP, σDOWN, σPLUS, σMINUS, σs, spinorrotation
using LinearAlgebra: I, det, tr

# Pauli-matrix algebra. These are tiny but catch any accidental sign flip,
# transpose swap, or rebinding of the σ constants.
@testset "Utils: Pauli matrices" begin
    @testset "definitions and aliases" begin
        @test σX === σ1
        @test σY === σ2
        @test σZ === σ3
        @test σ0 == [1.0 0.0; 0.0 1.0]
    end

    @testset "involution: σᵢ² = I" begin
        for σ in (σ1, σ2, σ3)
            @test σ * σ ≈ I
        end
    end

    @testset "tracelessness: tr σᵢ = 0" begin
        for σ in (σ1, σ2, σ3)
            @test isapprox(tr(σ), 0; atol=eps())
        end
    end

    @testset "anticommutation: {σᵢ, σⱼ} = 2δᵢⱼ I" begin
        for (i, σi) in enumerate((σ1, σ2, σ3)),
            (j, σj) in enumerate((σ1, σ2, σ3))
            @test σi * σj + σj * σi ≈ 2 * (i == j) * Matrix{ComplexF64}(I, 2, 2)
        end
    end

    @testset "commutation: [σ₁, σ₂] = 2i σ₃ (and cyclic)" begin
        @test σ1 * σ2 - σ2 * σ1 ≈ 2im * σ3
        @test σ2 * σ3 - σ3 * σ2 ≈ 2im * σ1
        @test σ3 * σ1 - σ1 * σ3 ≈ 2im * σ2
    end

    @testset "spin projectors: σ↑ + σ↓ = I and orthogonality" begin
        @test σUP + σDOWN ≈ σ0
        @test σUP * σDOWN ≈ zeros(2, 2)
        @test σUP * σUP ≈ σUP                 # idempotent
        @test σDOWN * σDOWN ≈ σDOWN
    end
end

# spinorrotation must be unitary for any θ, n.
@testset "Utils: spinorrotation is unitary" begin
    for θ in (0.0, π/4, π/2, π, 1.234), n in ([0,0,1], [1,0,0], [0,1,0], [1,1,1])
        U = spinorrotation(θ, Float64.(n))
        @test U * U' ≈ Matrix{ComplexF64}(I, 2, 2)
        @test isapprox(abs(det(U)), 1; atol=1e-12)
    end
end
