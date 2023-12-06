using Distributed
nprocs() < 2 && addprocs(47)

using LinearAlgebra
BLAS.set_num_threads(1)

include("system.jl")

using FileIO, DelimitedFiles

println("Building system...")
lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=true, spinhalf=true)

v = Operators.getcappedyukawa(lat; format=:dense, cellrange=3, spin=true, a=5.0, U=3.3) # interaction potential
ρ_init = Meanfield.initialguess(v, :random, :nonlocal; lat=lat) # initial guess

# Chemical potentials to sweep over
μs = LinRange(-0.09,-0.16,20)
fillings = Spectrum.filling(H, μs; nk=100)

mkpath("output_mf_sweep")
writedlm("output_mf_sweep/chemicalpotentials.out", μs)
writedlm("output_mf_sweep/fillings.out", fillings)

# Do mean-field calculations
println("Starting calculation...")
for (i_, filling) in enumerate(fillings)

    println("  doing #$i_ out of "*string(length(fillings)))

    Operators.setfilling!(H, filling; nk=120^2)

    ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock( # run the calculation
        H, v, ρ_init, filling; klin=120, iterations=1000, tol=1e-5,# p_norm=Inf,
        T=0.002, β=0.35,  show_trace=true, clear_trace=true, multimode=:distributed
    )

    # global ρ_init = ρ_sol # should accelarate convergence, but might also cause problems near phase transitions


    FileIO.save("output_mf_sweep/meanfield$i_.jld2", 
        Dict("ρ_sol"=>ρ_sol, "ϵ_GS"=>ϵ_GS, "HMF"=>HMF, "converged"=>converged, "error"=>error)
    )
end

println("Done.")
