using Distributed
nprocs() < 2 && addprocs(7)

using LinearAlgebra
BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "..", "system.jl"))
include(joinpath(@__DIR__, "params.jl"))

using FileIO, DelimitedFiles

println("Building system...")
lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=true, spinhalf=true)

v = Operators.getcappedyukawa(lat; format=:dense, cellrange=3, spin=true,
                                a=sp.a, U=sp.U)
ρ_init = Meanfield.initialguess(v, :random, :nonlocal; lat=lat)

# Chemical potentials to sweep over.
μs = LinRange(-0.09,-0.16,20)
fillings = Spectrum.filling(H, μs; nk=100)

mkpath(OUTDIR)
writedlm(joinpath(OUTDIR, "chemicalpotentials.out"), μs)
writedlm(joinpath(OUTDIR, "fillings.out"), fillings)
FileIO.save(joinpath(OUTDIR, "interaction_params.jld2"),
            Dict("a"=>sp.a, "U"=>sp.U))

iter = sp.direction === :reverse ? reverse(collect(enumerate(fillings))) :
                                    collect(enumerate(fillings))

println("Starting calculation...")
for (i_, filling) in iter
    println("  doing #$i_ out of $(length(fillings))")

    Operators.setfilling!(H, filling; nk=sp.nk^2)

    ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock(
        H, v, ρ_init, filling; klin=sp.klin, iterations=sp.iterations, tol=1e-5,
        T=0.002, β=sp.β, show_trace=true, clear_trace=true, multimode=:distributed
    )

    # Warm-start the next μ from the converged ρ. Speeds up convergence
    # but risks landing in the wrong basin near a phase transition; turn
    # off via `update_init=false` for sharper transition resolution.
    if sp.update_init
        global ρ_init = ρ_sol
    end

    FileIO.save(joinpath(OUTDIR, "meanfield$i_.jld2"),
        Dict("ρ_sol"=>ρ_sol, "ϵ_GS"=>ϵ_GS, "HMF"=>HMF,
             "converged"=>converged, "error"=>error)
    )
end

println("Done.")
