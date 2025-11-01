
"""
    fermisurfacedensity_fromdata(bands, energy; broadening=:auto)

Compute a broadened Fermi‑surface “density” from a band‑energy matrix `bands`
of size `(nbands, Nk)` evaluated on a k‑grid. `energy` can be a scalar or a
vector of Fermi energies. The result has size `(length(energies), Nk)` and sums
contributions from all bands using a Lorentzian of width `broadening`.

If `broadening == :auto` (default), a heuristic width based on the band
dispersion is used.
"""
fermisurfacedensity_fromdata(M, energy::Number, args...; kwargs...) = fermisurfacedensity_fromdata(M, [energy], args...; kwargs...)

"""
    getbroadening(broadening, bands)

Internal helper that resolves `broadening` to a numeric value. For `:auto`,
uses the mean band standard deviation divided by √Nk.
"""
function getbroadening(broadening, bands)
   if !isa(broadening, Number) || broadening==:auto
        broadening = Statistics.mean([Statistics.std(band) for band in eachrow(bands)])/sqrt(size(bands,2))
        println("Automatic broadening: $broadening")
   end

   broadening
end

"""
    fermisurfacedensity_fromdata(bands, energies; broadening=:auto)

Broadened Fermi‑surface density without weights/observables. See
`fermisurfacedensity_fromdata(bands, energy; ...)` for details.
"""
function fermisurfacedensity_fromdata(bands::AbstractMatrix, energies::AbstractVector; broadening=:auto)

    broadening::Float64 = getbroadening(broadening, bands)

    density = Array{Float64}(undef, (length(energies), size(bands,2)))

    for (j_, es) = enumerate(eachcol(bands))
        density[:, j_] .=  sum((1/pi) .* broadening^2 ./ (broadening^2 .+ (es.-energies').^2); dims=1)[:]
    end

    density
end

"""
    fermisurfacedensity_fromdata(bands, energies, obs; broadening=:auto)

Weighted variant where `obs` provides an observable per band and k‑point of the
same shape as `bands`. Each Lorentzian contribution is multiplied by the
corresponding weight.
"""
function fermisurfacedensity_fromdata(bands::AbstractMatrix, energies::AbstractVector, obs::AbstractMatrix; broadening=:auto)

    broadening::Float64 = getbroadening(broadening, bands)

    density = Array{Float64}(undef, (length(energies), size(bands,2)))

    for (j_, es) = enumerate(eachcol(bands))
        density[:, j_] .=  sum((1/pi) .* obs[:,j_] .* broadening^2 ./ (broadening^2 .+ (es.-energies').^2); dims=1)[:]
    end

    density
end

"""
    fermisurfacedensity(H, energies; broadening=:auto, lat=nothing, num_points=15, kwargs...)

High‑level convenience that computes band energies on a regular k‑grid and then
returns `(kgrid, density)` as produced by `fermisurfacedensity_fromdata`. If
`lat` is provided, the grid is folded into the first Brillouin zone and mapped
to Cartesian reciprocal space via `getB(lat)`.
"""
function fermisurfacedensity(H, args...; broadening=:auto, lat=nothing, num_points=15, kwargs...)

    kgrid = regulargrid(; nk=num_points^2)
    bands = bandmatrix(H, kgrid; kwargs...)[1]

    density = fermisurfacedensity_fromdata(bands, args...; broadening=broadening)

    if lat != nothing
        foldBZ!(kgrid, lat)
        kgrid = getB(lat) * kgrid
    end

    kgrid, density
end
