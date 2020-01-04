
using Distributed
using SharedArrays


@fastmath function ρ_k!(ρ0::T1, ϵs::T2, U::T3, μ::Float64; T::Float64=0.01, φk::ComplexF64=1.0+0.0im) where {T1<:AbstractMatrix{ComplexF64}, T2<:AbstractVector{ComplexF64}, T3<:AbstractMatrix{ComplexF64}}
    # ϵs, U = spectrum_k
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        ρ0[:] .+= (fermidirac(ϵ-μ; T=T) .* (ψ * ψ') .* φk)[:]
    end

    nothing
end

function ρ!(ρ0::T1, spectrum::Function, ks::T2, μ::Float64=0.0; T::Float64=0.01) where {T1<:AbstractMatrix{ComplexF64}, T2<:AbstractMatrix{Float64}}
    ρ0[:] .= zero(ρ0)[:]
    L = size(ks)[2]

    ρ0[:] = @distributed (+) for j=1:L
        ρtmp = zero(ρ0)    ## <-- it annoys me that I don't know how to get around this allocation

        ϵs, U = spectrum(ks[:,j])
        ρ_k!(ρtmp, ϵs, U, μ; T=T) #; format=format)

        ρtmp
    end

    ρ0[:] .= ρ0[:] ./ L

    nothing
end

function ρ(spectrum::Function, d::Int, ks::T, μ::Float64; kwargs...) where {T<:AbstractMatrix{Float64}}
    ρ0 = zeros(ComplexF64, d, d)
    ρ!(ρ0, spectrum, ks, μ; kwargs...)

    ρ0
end

################################################################################
################################################################################
################################################################################
################################################################################

using ProgressMeter

function ρ_L!(ρs::Dict{Vector{Int},T1}, spectrum::Function, ks::T2, μ::Float64=0.0; T::Float64=0.01) where {T1<:AbstractMatrix{ComplexF64}, T2<:AbstractMatrix{Float64}}
    L = size(ks,2)

    energies0_k = convert(SharedArray, zeros(Float64, L))

    for (δL,ρ0)=ρs
        ρs[δL][:] .= 0.0 #convert(SharedArray, zero(ρ0))[:]
    end

    @sync @distributed for i_=1:L

        k = ks[:,i_]
        energies_k, U_k = spectrum(k) #@time

        for δL=keys(ρs)

            for (ϵ, ψ) in zip(energies_k, eachcol(U_k)) # this loop used to be seperate in ρ_k!(...)
                ρs[δL][:] += (fermidirac(ϵ-μ; T=T) .* transpose(ψ * ψ') .* BlochPhase(-k, δL))[:]
            end

        end

        energies0_k[i_] = groundstate_sumk(energies_k, μ)
    end

    for δL = keys(ρs)
        ρs[δL][:] ./= L
    end

    sum(energies0_k)/L # return the groundstate energy
end
