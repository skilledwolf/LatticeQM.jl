module meanfieldtest
using Test
using LatticeQM
using LinearAlgebra: dot, norm
include("analyticbands.jl")
using .analyticbands: grapheneconductionband
using Plots

function get_graphene_HMF(U=4.0; filling=0.5, init=:antiferro, T=0.01, β=0.20)
    lat = Geometries.honeycomb()
    hops = Operators.graphene(lat; mode=:spinhalf)

    v = Operators.gethubbard(lat; mode=:σx, a=0.5, U=U)
    
    ρ_init = Meanfield.initialguess(v, init; lat=lat)

    ρ, _, HMF, converged, _ = Meanfield.solvehartreefock(
    hops, v, ρ_init, filling; klin=30, iterations=800, tol=1e-5,
    T=T, β=β,  show_trace=false, clear_trace=false)

    return lat, ρ, HMF, converged
end

function meanfieldsanitycheck()
    
    lat, _, HMF, converged = get_graphene_HMF(0.0; T=0.0, β=0.70)
    
    if !converged
        return false
    end

    ks = kpath(lat; num_points=200)
    nk = size(ks)[2]
    bands_mf = getbands(HMF.h, ks)
    bands_mf.bands .-= HMF.μ
    
    conductionband = map(grapheneconductionband, eachcol(ks.points))
    valenceband = - conductionband
    
    for band in eachrow(bands_mf.bands)
        maeavg1 = sum(broadcast(abs, band - conductionband)) / nk
        maeavg2 = sum(broadcast(abs, band - valenceband)) / nk

        # Small but significant differences from the analytical result should are expected
        # due to floating point impercision during calculation of the density matrix.
        if !(maeavg1 < 6e-4) && !(maeavg2 < 6e-4) 
            return false
        end
    end

    return true
end

function get_gap_at_U(U=4.0; filling=0.5, init=:antiferro, T=0.01, β=0.20)
    
    lat, _, HMF, _ = get_graphene_HMF(U; filling=filling, init=init, T=T, β=β)
    ks  = kpath(lat; num_points=200)
    gap = Spectrum.bandgap_filling(HMF.h, ks, filling)

    return gap
end

function hubbardmodelcriticalpoint()
    Us = range(0; stop=4.0, length=10)
    _, minidx = findmin(U -> abs(U - 2.23), Us)
    gaps = map(U -> get_gap_at_U(U; init=:antiferro), Us)

    dgaps = []
    for i in 2:length(gaps)
        push!(dgaps, gaps[i] - gaps[i - 1])
    end

    ddgaps = []
    for i in 2:length(dgaps)
        push!(ddgaps, dgaps[i] - dgaps[i - 1])
    end

    maxdd = max(ddgaps...)

    # Critical value should be defined as the U value at which the second derivative
    # of the band gap opening finds it's maximum.
    return ddgaps[minidx - 2] == maxdd || ddgaps[minidx - 1] == maxdd
end

function hubbardmodelantiferro()
    lat, ρ, _, _ = get_graphene_HMF()
    suba = getoperator(lat, "sublatticeAspin")
    subb = getoperator(lat, "sublatticeBspin")

    subamag_ops = [suba*s for s=getoperator(lat, ["sx", "sy", "sz"])]
    subbmag_ops = [subb*s for s=getoperator(lat, ["sx", "sy", "sz"])]
    subamag = real([Operators.trace(ρ, mag) for mag in subamag_ops])
    subbmag = real([Operators.trace(ρ, mag) for mag in subbmag_ops])
    
    return isapprox(dot(subamag, subbmag) / norm(subamag) / norm(subbmag), -1.0; atol=1e-12)
end

@test meanfieldsanitycheck()
@test hubbardmodelcriticalpoint()
@test hubbardmodelantiferro()
end