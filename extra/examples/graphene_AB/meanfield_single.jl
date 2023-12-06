include("system.jl")

lat, H, ops = buildsystem(; tz=0.3, V=0.3, triple=true, spinhalf=true)

μ = -0.11 #-0.14 #-0.105 # energies[3]

filling = Spectrum.filling(H, μ; nk=100)
Operators.setfilling!(H, filling; nk=100^2)

v = Operators.getcappedyukawa(lat; format=:dense, cellrange=3, spin=true, a=4.0, U=3.0) # interaction potential
ρ_init = Meanfield.initialguess(v, :random, :nonlocal; lat=lat) # initial guess
# ρ_init = Meanfield.initialguess(v, :ferro; lat=lat) # initial guess

ρ_sol, ϵ_GS, HMF, converged, error = Meanfield.solvehartreefock( # run the calculation
    H, v, ρ_init, filling; klin=50, iterations=350, tol=1e-5,# p_norm=Inf,
    T=0.002, β=0.45,  show_trace=true, clear_trace=true
)


using FileIO
mkpath("ouput_mf")
save("meanfield1.jld2", "ρ_sol", ρ_sol, "ϵ_GS", ϵ_GS, "HMF", "converged", converged, "error", error)