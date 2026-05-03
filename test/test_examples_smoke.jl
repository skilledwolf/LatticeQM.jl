using Test
using LatticeQM

# Smoke tests for `extra/examples/`: include each fast example end-to-end
# and assert it runs to completion. Catches the kind of API rot we hit
# when `HartreeFock(::SparseHops)` started erroring (the bug was real;
# nothing in the test suite exercised the sparse + HartreeFock path before
# this file existed).
#
# Gated on `LATTICEQM_FULL_EXAMPLES=1` so the default `Pkg.test()` stays
# fast and avoids a Plots dependency. Most examples `using Plots`, so we
# also skip if Plots is not installed in the active environment.

const FULL_EXAMPLES = get(ENV, "LATTICEQM_FULL_EXAMPLES", "") == "1"
const EXAMPLES_DIR  = abspath(joinpath(@__DIR__, "..", "extra", "examples"))

# Examples that exercise the public API in a few seconds and don't open
# windows or require external data. Slow examples (SCF sweeps, Hofstadter,
# distributed-only ones) are deliberately excluded; they're tested
# implicitly by `test_meanfield.jl` etc.
const FAST_EXAMPLES = [
    "graphene/dos.jl",
    "graphene/valleybands.jl",
    "graphene/ribbons.jl",
    "graphene/ribbon_haldane.jl",
    "graphene/opticalconductivity.jl",
    "graphene/hubbardmeanfield_example2.jl",
    "twistedgraphene/valleybands.jl",
    "SSHxSSH/wilsonloop.jl",
]

# Each example is run in a fresh `Module` so its `const` / top-level
# definitions don't leak between tests, and inside `mktempdir() / cd()`
# so output files land somewhere disposable.
function run_example(rel_path::String)
    src = joinpath(EXAMPLES_DIR, rel_path)
    @assert isfile(src) "missing example: $rel_path"
    mktempdir() do tmp
        prev = pwd()
        cd(tmp)
        try
            mod = Module(:ExampleSmoke)
            Base.invokelatest(() -> Base.include(mod, src))
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

    # Most examples `using Plots`. Probe once up front; bail with a clear
    # message rather than erroring on every example.
    plots_ok = try
        @eval Main using Plots
        true
    catch err
        @warn "Plots not available — examples smoke tests will be skipped" err
        false
    end
    plots_ok || return

    ENV["GKSwstype"] = "nul"  # GR headless

    for rel in FAST_EXAMPLES
        @testset "$rel" begin
            try
                run_example(rel)
                @test true
            catch err
                @error "Example failed" rel err
                @test false
            end
        end
    end
end
