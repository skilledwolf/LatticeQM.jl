using Plots, LinearAlgebra
using LatticeQM
using LatticeQM.Operators: graphene, addzeeman!, magnetization, density
using LatticeQM.Meanfield

roundreal(x; digits=7) = round.(real.(x); digits=digits)

function get_gap_at_U(U=4.0; filling=0.5, init=:antiferro, T=0.01, β=0.20, show_trace=true, show_bands=false, clear_trace=false, reportmagnetization=false)

    lat = Geometries2D.honeycomb()
    sx, sy, sz, sublA, sublB = getoperator(lat, ["sx", "sy", "sz", "sublatticeAspin", "sublatticeBspin"])
    hops = graphene(lat; mode=:spinhalf)
#     addzeeman!(hops, lat, 0.0001)

    v = gethubbard(lat; mode=:σx, a=0.5, U=U) # interaction potential
    ρ_init = initialguess(v, init; lat=lat) # initial guess

    hf = hartreefock(hops, v)

    ρ_sol, ϵ_GS, HMF, converged, error = solveselfconsistent(
        hf, ρ_init, filling; klin=30, iterations=500, tol=1e-8,# p_norm=Inf,
        T=T, β=β,  show_trace=show_trace, clear_trace=clear_trace
    )

    ks  = kpath(lat; num_points=200)
    gap = Spectrum.bandgap_filling(HMF.h, ks, filling)

    if show_bands
        bands = getbands(hops, ks, sz)
        p1 = plot(bands; size=(400,300))

        bandsMF = getbands(HMF.h, ks, sz)
        bandsMF.bands .-= HMF.μ

        p2 = plot(bandsMF; size=(400,300))
        p = plot(p1,p2, size=(800,300)) #gui(...)
        display(p)
    end

    if reportmagnetization
        mA, mB = roundreal.(magnetization(ρ_sol, [sublA,sublB], lat))
        δM = mA - mB
        M = mA+mB
        Mabs = norm(M)
        δMabs = norm(δM)
        dens = density(ρ_sol)
        @info("Groundstate energy", ϵ_GS)
        @info("Magnetization", Mabs, δMabs, mA, mB, M, δM)
        @info("Density", dens)
    end

    [gap ϵ_GS]
end

get_gap_at_U(3.0; filling=0.5, init=:antiferro, T=0.01, β=0.15, show_bands=true, clear_trace=true, reportmagnetization=true)
mkpath("output"); savefig("output/bandsmf.pdf")

using ProgressMeter

Us = range(0; stop=4.0, length=10)
data = progress_map(Us) do U
    get_gap_at_U(U; init=:antiferro, show_trace=false)
end

data = vcat(data...)
gap = data[:,1]
energies = data[:,2]

scatter(Us, gap; ylim=(0,3), size=(300,200), xlabel="U", ylabel="band gap")
mkpath("output"); savefig("output/hubbardsweep_gap.pdf")


