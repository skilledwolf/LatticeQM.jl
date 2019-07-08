using Distributed

function dos_dense(hamiltonian::Function, ks::AbstractMatrix{Float64}, energies::AbstractVector{Float64}; Γ::Float64)

    L = size(ks)[2]
    dos = zero(energies)

    for j=1:L
        for ϵ in ϵs_dense(hamiltonian)(ks[:,j])
            dos .+= imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) )
        end
    end

    dos / L / π
end

function dos_dense_parallel(hamiltonian::Function, ks::AbstractMatrix{Float64}, energies::AbstractVector{Float64}; Γ::Float64)

    L = size(ks)[2]

    dos = @distributed (+) for j=1:L # over ks
        tmpdos = zero(energies)

        for ϵ in ϵs_dense(hamiltonian)(ks[:,j]) # over bands at k
            tmpdos .+= imag.( 1.0./(energies .- 1.0im .* Γ .- ϵ) )
        end

        tmpdos / L / π
    end

    dos
end

function dos_dense_parallel(hamiltonian::Function, ks::AbstractMatrix{Float64}, ωmin::Float64, ωmax::Float64, num::Int)

    energies = collect(range(ωmin, length=num, stop=ωmax))
    Γ = 0.1 * abs(ωmax-ωmin)/(num-1)

    energies, dos_dense_parallel(hamiltonian, ks, energies; Γ=Γ)
end
