import ..Utils: regulargrid

function bandgap_filling(H, filling::Real; klin=30, kwargs...)
    kgrid = regulargrid(;nk=klin^2)
    bandgap_filling(H, kgrid, filling; kwargs...)
end

function bandgap_filling(H, ks, filling::Real; multimode=:distributed, kwargs...)
    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(H, ks; multimode=multimode) # dense diagonalization (default)!
    μ = chemicalpotential(bands, filling; kwargs...)

    bandgap_energy(bands, μ)
end

function bandgap(H, μ::Real=0.0; klin=10, multimode=:distributed) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    kgrid = regulargrid(;nk=klin^2)
    bandgap_energy(bandmatrix(H, kgrid; multimode=multimode), μ)
end

function bandgap_energy(H, ks, μ::Real=0.0, multimode=:distributed) # Note: dense diagonalization!
    # Calculate the gap around in which the Fermi level lies
    bandgap_energy(bandmatrix(H, ks; multimode=multimode), μ)
end

function bandgap_energy(bands::AbstractMatrix, μ::Real)

    ϵlower = maximum(bands[bands.<= μ])
    ϵupper = minimum(bands[bands.>= μ])
    ϵgap = ϵupper - ϵlower

    ϵgap
end