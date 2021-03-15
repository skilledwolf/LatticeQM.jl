import ..Utils: regulargrid

function bandgap_filling(H, filling::Real; klin=30, kwargs...)
    kgrid = regulargrid(;nk=klin^2)
    bandgap_filling(H, kgrid, filling; kwargs...)
end

function bandgap_filling(H, ks, filling::Real; kwargs...)
    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(H, points(ks)) # dense diagonalization (default)!
    μ = chemicalpotential(bands, filling; kwargs...)

    bandgap_energy(bands, μ)
end

function bandgap(H, μ::Real=0.0; klin=10) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    kgrid = regulargrid(;nk=klin^2)
    bandgap_energy(bandmatrix(H, kgrid), μ)
end

function bandgap_energy(H, ks, μ::Real=0.0) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    bandgap_energy(bandmatrix(H, points(ks)), μ)
end

function bandgap_energy(bands::AbstractMatrix, μ::Real)

    ϵlower = maximum(bands[bands.<= μ])
    ϵupper = minimum(bands[bands.>= μ])
    ϵgap = ϵupper - ϵlower

    ϵgap
end