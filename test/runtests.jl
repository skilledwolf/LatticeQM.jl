using Test

# Top-level test driver. Each include defines its own @testset(s); this file
# wraps them in a single named outer set so failures attribute cleanly.
#
# The default `Pkg.test()` runs the fast unit + regression suite (~15 s).
# Set LATTICEQM_FULL_TESTS=1 to also run the slow Hartree-Fock SCF regression
# (~5 min) — appropriate for nightly/release CI but not every push.
@testset "LatticeQM" begin
    @testset "Utils"               include("test_utils.jl")
    @testset "Structure"           include("test_structure.jl")
    @testset "TightBinding"        include("test_tightbinding.jl")
    @testset "Eigen"               include("test_eigen.jl")
    @testset "Operators"           include("test_operators.jl")
    @testset "Spectrum/graphene"   include("test_spectrum_graphene.jl")
    @testset "Spectrum/square"     include("test_spectrum_square.jl")
    @testset "Spectrum/triangular" include("test_spectrum_triangular.jl")
    @testset "Superconductivity"   include("test_superconductivity.jl")
    @testset "Floquet"             include("test_floquet.jl")
    @testset "Haldane"             include("test_haldane.jl")
    @testset "TBG"                 include("test_tbg.jl")
    @testset "Meanfield"           include("test_meanfield.jl")
    @testset "Parallelism"         include("test_parallelism.jl")
    @testset "Examples (smoke)"    include("test_examples_smoke.jl")
end
