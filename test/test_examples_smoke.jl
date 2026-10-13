using Test
using LatticeQM

# Smoke tests for `extra/examples/`: include each fast example end-to-end
# and assert it runs to completion. Catches the kind of API rot we hit
# when `HartreeFock(::SparseHops)` started erroring (the bug was real;
# nothing in the test suite exercised the sparse + HartreeFock path before
# this file existed).
#
# Gated on `LATTICEQM_FULL_EXAMPLES=1` so the default `Pkg.test()` stays
# fast and avoids a Plots dependency. Explicitly enabling these tests requires
# Plots in the active environment; a missing dependency must fail the run.

const FULL_EXAMPLES = get(ENV, "LATTICEQM_FULL_EXAMPLES", "") == "1"
const EXAMPLES_DIR  = abspath(joinpath(@__DIR__, "..", "extra", "examples"))

# Examples that exercise the public API without opening windows or requiring
# external data. Some retain large grids to check realistic workloads.
# Longer examples (SCF sweeps, Hofstadter,
# distributed-only ones) are deliberately excluded; they're tested
# implicitly by `test_meanfield.jl` etc.
const FAST_EXAMPLES = [
    "graphene/dos.jl" => ("dos.pdf",),
    "graphene/valleybands.jl" => ("bands.pdf",),
    "graphene/ribbons.jl" => ("ribbons_geometry.pdf", "ribbons_bands.pdf"),
    "graphene/ribbon_haldane.jl" => ("ribbon_haldane_bands.pdf",),
    "graphene/opticalconductivity.jl" => ("bands_oc.pdf", "opticalconductivity.pdf"),
    "graphene/hubbardmeanfield_example2.jl" => ("hubbardmeanfield_example2.pdf",),
    "twistedgraphene/valleybands.jl" => ("valleybands.pdf",),
    "SSHxSSH/wilsonloop.jl" => ("pol1.pdf", "pol2.pdf"),
]

# Each example is run in a fresh `Module` so its `const` / top-level
# definitions don't leak between tests, and inside `mktempdir() / cd()`
# so output files land somewhere disposable.
function run_example(rel_path::String, expected_pdfs)
    src = joinpath(EXAMPLES_DIR, rel_path)
    @assert isfile(src) "missing example: $rel_path"
    mktempdir() do tmp
        prev = pwd()
        cd(tmp)
        try
            mod = Module(:ExampleSmoke)
            Base.invokelatest(() -> Base.include(mod, src))
            for filename in expected_pdfs
                path = joinpath(tmp, "output", filename)
                @test isfile(path) && filesize(path) > 0
            end
        finally
            cd(prev)
        end
    end
    nothing
end

@testset "examples smoke" begin
    if !FULL_EXAMPLES
        @info "Skipping examples smoke tests (set LATTICEQM_FULL_EXAMPLES=1 to run)"
        return
    end

    # Configure GR before loading Plots, and restore the caller's environment
    # even if a dependency or example fails.
    withenv("GKSwstype" => "nul") do
        @eval Main using Plots

        for (rel, expected_pdfs) in FAST_EXAMPLES
            @testset "$rel" begin
                @info "Running example smoke test" example=rel
                started = time()
                try
                    run_example(rel, expected_pdfs)
                finally
                    @info "Finished example smoke test" example=rel seconds=round(time() - started; digits=2)
                end
            end
        end
    end
end
