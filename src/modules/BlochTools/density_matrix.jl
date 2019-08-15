
using Distributed
using SharedArrays


function ρ_k!(ρ0::AbstractMatrix{ComplexF64}, spectrum_k, μ::Float64; T::Float64=0.01, φk::ComplexF64=1.0+0.0im)
    ϵs, U = spectrum_k
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        ρ0[:] .+= (fermidirac(ϵ-μ; T=T) .* (ψ * ψ') .* φk)[:]
    end

    nothing
end

function ρ!(ρ0::AbstractMatrix{ComplexF64}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)
    ρ0[:] .= zero(ρ0)[:]
    L = size(ks)[2]

    ρ0[:] = @distributed (+) for j=1:L
        ρtmp = zero(ρ0)    ## <-- it annoys me that I don't know how to get around this allocation
        ρ_k!(ρtmp, spectrum(ks[:,j]), μ) #; format=format)

        ρtmp
    end

    ρ0[:] .= ρ0[:] ./ L

    nothing
end

function ρ(spectrum::Function, d::Int, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    ρ0 = zeros(ComplexF64, d, d)
    ρ!(ρ0, spectrum, ks, μ; kwargs...)

    ρ0
end

################################################################################
################################################################################
################################################################################
################################################################################

using ProgressMeter

function ρ_L!(ρs::Dict{Vector{Int},AbstractMatrix{ComplexF64}}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; T::Float64=0.01)
    L = size(ks,2)

    energies0_k = convert(SharedArray, zeros(Float64, L))

    for (δL,ρ0)=ρs
        ρs[δL][:] .= convert(SharedArray, zero(ρ0))[:]
    end

    ## EXPERIMENTAL USAGE OF THE PROGRESS BAR
    p = Progress(L, 0.1, "Computing density matrix...")
    channel = RemoteChannel(()->Channel{Bool}(L), 1)

    @sync begin
        # this task prints the progress bar
        @async while take!(channel)
            next!(p)
        end

        # this task does the computation
        @async begin
            @distributed for i_=1:L
            # for i_=1:L
                k = ks[:,i_]
                spectrum_k = spectrum(k)

                for δL=keys(ρs)
                    ρ_k!(ρs[δL], spectrum_k, μ; φk=BlochPhase(-k, δL), T=T)
                end

                energies0_k[i_] = groundstate_sumk(spectrum_k[1], μ)
            end

            put!(channel, false) # this tells the printing task to finish
        end
    end
    ## END OF EXPERIMENTAL USE OF THE PROGRESSBAR

    for δL = keys(ρs)
        ρs[δL] ./= L
    end

    sum(energies0_k)/L # return the groundstate energy
end
