using Test
using LatticeQM
using LinearAlgebra: diag, tr, Diagonal, ishermitian

# J-augmented Hartree-Fock (Meanfield.ExchangeKernel): analytic limits,
# variational consistency, hermiticity, and ferromagnetic selection.

# ---------------------------------------------------------------- helpers
# Two sites x two flavors per cell, layout [s1f1, s2f1, s1f2, s2f2]
# (flavor-outer), on a 1D-like lattice with displacements [0,0] and [±1,0].
const SITEOF = [1, 2, 1, 2]
const FLAVOROF = [1, 1, 2, 2]

function make_jkernel(J1::Float64, Jnn::Float64)
    # intra-cell exchange between sites 1-2 (J1), inter-cell site1-site1 (Jnn)
    j0 = [0.0 J1; J1 0.0]
    jp = [Jnn 0.0; 0.0 0.0]
    Hops([0, 0] => complex(j0), [1, 0] => complex(jp),
         [-1, 0] => complex(jp'))
end

function random_hermitian_rho(D)
    ρ0 = rand(ComplexF64, D, D)
    ρ0 = (ρ0 + ρ0') / 2
    ρ1 = rand(ComplexF64, D, D)
    Hops([0, 0] => ρ0, [1, 0] => ρ1, [-1, 0] => collect(ρ1'))
end

@testset "ExchangeKernel construction" begin
    ek = ExchangeKernel(make_jkernel(1.0, 0.3), SITEOF, FLAVOROF)
    @test size(ek.sf2idx) == (2, 2)
    @test ek.sf2idx[1, 1] == 1 && ek.sf2idx[2, 2] == 4
    # self-exchange on the diagonal of the zero key must be rejected
    jbad = Hops([0, 0] => complex([0.5 0.0; 0.0 0.0]))
    @test_throws ArgumentError ExchangeKernel(jbad, SITEOF, FLAVOROF)
    # mismatched site count must be rejected
    @test_throws ArgumentError ExchangeKernel(make_jkernel(1.0, 0.0),
                                              [1, 1, 1, 1], FLAVOROF)
end

@testset "analytic dimer limit" begin
    # J only between sites 1 and 2 in the same cell; diagonal ρ.
    J1 = 0.7
    h = Hops([0, 0] => complex(zeros(4, 4)))
    v = Hops([0, 0] => complex(0.0 * ones(4, 4)))
    ek = ExchangeKernel(make_jkernel(J1, 0.0), SITEOF, FLAVOROF)
    hf = HartreeFock(h, v; exchange=ek)

    # fully polarized: both sites filled in flavor 1 -> eJ = -J1
    ρpol = Hops([0, 0] => complex(Matrix(Diagonal([1.0, 1.0, 0.0, 0.0]))))
    hf(ρpol)
    @test isapprox(hf.ϵJ, -J1; atol=1e-12)

    # unpolarized: n = 1/2 per flavor -> eJ = -J1/2
    ρunp = Hops([0, 0] => complex(Matrix(Diagonal([0.5, 0.5, 0.5, 0.5]))))
    hf(ρunp)
    @test isapprox(hf.ϵJ, -J1 / 2; atol=1e-12)
end

@testset "variational identity and hermiticity" begin
    # Tr[Σ_J ρ] = 2 ϵJ for generic complex Hermitian ρ, including bond keys.
    h = Hops([0, 0] => complex(zeros(4, 4)))
    v = Hops([0, 0] => complex(zeros(4, 4)))
    ek = ExchangeKernel(make_jkernel(0.8, 0.35), SITEOF, FLAVOROF)
    hf0 = HartreeFock(h, v)              # reference without exchange
    hfJ = HartreeFock(h, v; exchange=ek)

    ρ = random_hermitian_rho(4)
    hf0(ρ)
    hMF0 = deepcopy(hf0.hMF)
    hfJ(ρ)
    ΣJ = hfJ.hMF - hMF0                               # isolate the J fields
    @test TightBinding.ishermitian(ΣJ)
    @test isapprox(real(Operators.trace(ΣJ, ρ)), 2 * hfJ.ϵJ;
                   atol=1e-10, rtol=1e-10)
end

@testset "SCF: exchange selects the polarized state" begin
    # Two-site cell, two flavors, weak hopping between the sites, moderate
    # flavor-blind repulsion, half filling of one flavor sector (quarter
    # filling total). Without J the polarized and unpolarized solutions are
    # near-degenerate at this filling; with J the polarized state must win
    # by ~J per cell.
    t = 0.1
    U = 1.0
    J1 = 0.4
    h1 = complex([0.0 t; t 0.0])
    h = Hops([0, 0] => kron([1.0 0; 0 1.0], h1))   # flavor-outer layout
    vmat = complex(U * ones(4, 4) - U * Diagonal(ones(4)))  # offsite-ish
    v = Hops([0, 0] => vmat)
    ek = ExchangeKernel(make_jkernel(J1, 0.0), SITEOF, FLAVOROF)

    occs_pol = [0.9, 0.9, 0.1, 0.1]
    occs_unp = [0.5, 0.5, 0.5, 0.5]
    Es = Dict{String,Float64}()
    for (label, occs, exch) in [("pol_J", occs_pol, ek),
                                ("unp_J", occs_unp, ek),
                                ("pol_0", occs_pol, nothing),
                                ("unp_0", occs_unp, nothing)]
        ρ0 = Hops([0, 0] => complex(Matrix(Diagonal(occs))))
        ρ, eGS, mf, conv, res = Meanfield.solvehartreefock(
            h, v, ρ0, 0.5; klin=8, iterations=600, tol=1e-9, T=0.02,
            β=0.5, show_trace=false, exchange=exch)
        @test conv
        Es[label] = eGS
        # total-energy identities with the ϵJ bookkeeping
        @test isapprox(eGS, mf.ϵkin + mf.ϵH + mf.ϵF + mf.ϵJ;
                       atol=1e-6, rtol=1e-8)
    end
    # exchange gain of polarization ~ J1 (both states feel U equally here)
    gain_J = Es["unp_J"] - Es["pol_J"]
    gain_0 = Es["unp_0"] - Es["pol_0"]
    @test gain_J - gain_0 > 0.3 * J1
end
