function bandgap_filling_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, filling::Float64; kwargs...)

    # Calculate the gap around in which the Fermi level lies
    bands = bandmatrix(hamiltonian, ks)
    μ = chemicalpotential(bands, ks, filling; kwargs...)

    bandgap_μ_dense(bands, μ)
end

function bandgap_μ_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64)

    # Calculate the gap around in which the Fermi level lies
    bandgap_μ_dense(bandmatrix(hamiltonian, ks), μ)
end

function bandgap_μ_dense(bands::AbstractMatrix, μ::Float64)

    ϵlower = maximum(bands[bands.<= μ])
    ϵupper = minimum(bands[bands.>= μ])
    ϵgap = ϵupper - ϵlower

    ϵgap
end