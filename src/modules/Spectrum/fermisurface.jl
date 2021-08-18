


fermisurfacedensity_fromdata(M, energy::Number, args...; kwargs...) = fermisurfacedensity_fromdata(M, [energy], args...; kwargs...)

function getbroadening(broadening, bands)
   if !isa(broadening, Number) || broadening==:auto
        broadening = Statistics.mean([Statistics.std(band) for band in eachrow(bands)])/sqrt(size(bands,2))
        println("Automatic broadening: $broadening")
   end

   broadening
end

function fermisurfacedensity_fromdata(bands::AbstractMatrix, energies::AbstractVector; broadening=:auto)

    broadening::Float64 = getbroadening(broadening, bands)

    density = Array{Float64}(undef, (length(energies), size(bands,2)))

    for (j_, es) = enumerate(eachcol(bands))
        density[:, j_] .=  sum((1/pi) .* broadening^2 ./ (broadening^2 .+ (es.-energies').^2); dims=1)[:]
    end

    density
end

function fermisurfacedensity_fromdata(bands::AbstractMatrix, energies::AbstractVector, obs::AbstractMatrix; broadening=:auto)

    broadening::Float64 = getbroadening(broadening, bands)

    density = Array{Float64}(undef, (length(energies), size(bands,2)))

    for (j_, es) = enumerate(eachcol(bands))
        density[:, j_] .=  sum((1/pi) .* obs[:,j_] .* broadening^2 ./ (broadening^2 .+ (es.-energies').^2); dims=1)[:]
    end

    density
end



function fermisurfacedensity(H, args...; broadening=:auto, lat=nothing, num_points=15, kwargs...)

    kgrid = regulargrid(; nk=num_points^2)
    bands = bandmatrix(H, kgrid; kwargs...)

    density = fermisurfacedensity_fromdata(bands, args...; broadening=broadening)

    if lat != nothing
        foldBZ!(kgrid, lat)
        kgrid = getB(lat) * kgrid
    end

    kgrid, density
end