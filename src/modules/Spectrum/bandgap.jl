using ..Utils: regulargrid

function bandgap_filling(H, filling::Float64; klin=30, kwargs...)
    kgrid = regulargrid(;nk=klin^2)
    bandgap_filling(H, kgrid, filling; kwargs...)
end

function bandgap_filling(H, ks, filling::Float64; kwargs...)
    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(H, points(ks)) # dense diagonalization (default)!
    μ = chemicalpotential(bands, filling; kwargs...)

    bandgap_energy(bands, μ)
end

function bandgap_energy(H, ks, μ::Float64) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    bandgap_energy(bandmatrix(H, points(ks)), μ)
end

function bandgap_energy(bands::AbstractMatrix, μ::Float64)

    ϵlower = maximum(bands[bands.<= μ])
    ϵupper = minimum(bands[bands.>= μ])
    ϵgap = ϵupper - ϵlower

    ϵgap
end