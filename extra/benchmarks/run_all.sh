#!/bin/bash
# Drive bench_one.jl across {case} × {mode} and write a CSV.
set -euo pipefail
cd "$(dirname "$0")/../.."

OUT="${1:-extra/benchmarks/results.csv}"

echo "case,mode,nthreads,nworkers,blas,nks,N,time_s,rss_mb,allocs_mb" > "$OUT"

run() {
    local case="$1" mode="$2" threads="${3:-1}" extra="${4:-}"
    local cmd
    if [[ "$mode" == "threaded" ]]; then
        cmd="JULIA_NUM_THREADS=$threads julia --project=. --color=no extra/benchmarks/bench_one.jl $case $mode"
    else
        cmd="julia --project=. --color=no extra/benchmarks/bench_one.jl $case $mode $extra"
    fi
    echo ">> $case / $mode (threads=$threads, extra=$extra)" >&2
    eval "$cmd" 2>/dev/null | tee -a "$OUT"
}

for case in dense_small tbg_n3 tbg_n5 tbg_n7 tbg_n11; do
    run "$case" serial
    run "$case" threaded 8
    run "$case" distributed 1 4
done
