
# Sweep parameter sets, replacing the previous `graphene_AB_sweep1..5,1_reverse/`
# folders. Each entry sets the (a, U) point for the capped-Yukawa interaction
# plus the SCF tuning knobs. Pick one via `julia <script>.jl <name>` or
# `LATTICEQM_PARAMS=<name>`. Default is `:sweep1`.
#
# Direction `:forward` walks fillings high→low (or whatever the sense of
# `LinRange(...)` gives); `:reverse` walks the same fillings backward — the
# old `graphene_AB_sweep1_reverse` variant. Different orderings can land in
# different SCF basins near phase boundaries.

const SWEEP_PARAMS = Dict(
    :sweep1   => (a=5.0, U=3.3, klin=50,  iterations=500,  β=0.45, nk=100, update_init=true,  direction=:forward),
    :sweep1r  => (a=5.0, U=3.3, klin=50,  iterations=500,  β=0.45, nk=100, update_init=true,  direction=:reverse),
    :sweep2   => (a=4.0, U=3.6, klin=50,  iterations=900,  β=0.45, nk=100, update_init=true,  direction=:forward),
    :sweep3   => (a=5.0, U=3.3, klin=75,  iterations=900,  β=0.35, nk=100, update_init=true,  direction=:forward),
    :sweep4   => (a=5.0, U=3.3, klin=120, iterations=1000, β=0.35, nk=120, update_init=false, direction=:forward),
    :sweep5   => (a=3.0, U=3.5, klin=75,  iterations=900,  β=0.35, nk=100, update_init=true,  direction=:forward),
)

const params_name = Symbol(get(ENV, "LATTICEQM_PARAMS",
                                isempty(ARGS) ? "sweep1" : ARGS[1]))
haskey(SWEEP_PARAMS, params_name) ||
    error("Unknown sweep param set: $params_name. Available: $(collect(keys(SWEEP_PARAMS)))")
const sp = SWEEP_PARAMS[params_name]
@info "Running sweep param set" params_name sp

const OUTDIR = "output_$(params_name)"
