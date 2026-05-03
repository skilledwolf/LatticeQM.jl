# Compare two benchmark CSV files produced by `bench_one.jl`.
#
# Usage:  julia --project=. extra/benchmarks/diff_csv.jl <baseline.csv> <new.csv> [--threshold=0.10]
#
# Reports per-(case, mode) wall-time and allocation deltas. Exits non-zero
# (1) if any case is slower by more than `threshold` (default 10%) — which
# makes this script a drop-in CI guard against perf regressions.
#
# CSV format (matches `bench_one.jl`):
#   case,mode,nthreads,nworkers,blas,nks,N,time_s,rss_mb,allocs_mb

using DelimitedFiles

function read_results(path::AbstractString)
    raw, _ = readdlm(path, ','; header=true)
    rows = Dict{Tuple{String,String}, Dict{Symbol, Any}}()
    for r in eachrow(raw)
        key = (string(r[1]), string(r[2]))
        rows[key] = Dict(
            :nthreads => Int(r[3]),
            :nworkers => Int(r[4]),
            :blas     => Int(r[5]),
            :nks      => Int(r[6]),
            :N        => Int(r[7]),
            :time_s   => parse(Float64, string(r[8])),
            :rss_mb   => parse(Float64, string(r[9])),
            :allocs_mb => parse(Float64, string(r[10])),
        )
    end
    rows
end

function fmt_pct(x::Real)
    sign = x >= 0 ? "+" : ""
    string(sign, round(100*x; digits=1), "%")
end

function diff_table(baseline::AbstractString, new::AbstractString; threshold::Float64=0.10)
    base = read_results(baseline)
    cur = read_results(new)
    keys_all = sort(collect(union(keys(base), keys(cur))))

    println("# benchmark diff:  $baseline  →  $new   (regression threshold = $(round(100*threshold; digits=0))%)")
    println()
    println(rpad("case", 14), rpad("mode", 13),
            rpad("base_t", 10), rpad("new_t", 10), rpad("Δt", 10),
            rpad("base_alloc", 12), rpad("new_alloc", 12), rpad("Δalloc", 10))
    println(repeat("-", 92))

    regressed = String[]
    for k in keys_all
        b = get(base, k, nothing)
        c = get(cur, k, nothing)
        if b === nothing
            println(rpad(k[1], 14), rpad(k[2], 13), "(new only)")
            continue
        elseif c === nothing
            println(rpad(k[1], 14), rpad(k[2], 13), "(missing)")
            continue
        end
        Δt    = (c[:time_s] - b[:time_s]) / b[:time_s]
        Δa    = (c[:allocs_mb] - b[:allocs_mb]) / max(b[:allocs_mb], 1e-3)
        marker = Δt > threshold ? "  ⚠" : ""
        println(rpad(k[1], 14), rpad(k[2], 13),
                rpad(string(b[:time_s]), 10), rpad(string(c[:time_s]), 10), rpad(fmt_pct(Δt), 10),
                rpad(string(b[:allocs_mb]), 12), rpad(string(c[:allocs_mb]), 12), rpad(fmt_pct(Δa), 10),
                marker)
        if Δt > threshold
            push!(regressed, "$(k[1]) / $(k[2]) (+$(round(100*Δt; digits=1))%)")
        end
    end

    println()
    if isempty(regressed)
        println("OK: no regressions above $(round(100*threshold; digits=0))%")
        return 0
    else
        println("REGRESSIONS:")
        for r in regressed; println("  - $r"); end
        return 1
    end
end

function _main(args)
    threshold = 0.10
    files = String[]
    for a in args
        if startswith(a, "--threshold=")
            threshold = parse(Float64, a[13:end])
        else
            push!(files, a)
        end
    end
    length(files) == 2 || error("Usage: diff_csv.jl <baseline.csv> <new.csv> [--threshold=FRAC]")
    return diff_table(files[1], files[2]; threshold=threshold)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_main(ARGS))
end
