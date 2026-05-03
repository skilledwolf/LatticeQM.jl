using ProgressMeter

import ..Utils: fermidirac
using ..TightBinding: Hops, AbstractHops, dim
import LatticeQM.Parallel

function green(H, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    d = dim(H, ks)
    G = zeros(ComplexF64, d, d)
    green!(G, H, ks, μ; kwargs...)
    G
end

green(ϵ::Number, ψ::AbstractVector, ω::Number=0.0) = transpose(ψ * ψ') / (ω-ϵ) #transpose(ψ * ψ') # (ψ * ψ')

function green!(G::AbstractMatrix, ϵs::AbstractVector, U::AbstractMatrix; φk::ComplexF64=1.0+0.0im, kwargs...)
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        G[:] .+= (green(ϵ, ψ; kwargs...) .* φk)[:]
    end
    G
end

function green!(G::AbstractHops, k::AbstractVector, ϵs::AbstractVector, U::AbstractMatrix; kwargs...)
    for δL=keys(G)
        for (ϵ, ψ) in zip(ϵs, eachcol(U))
            G[δL][:] .+= (green(ϵ, ψ; kwargs...) .* fourierphase(-k, δL))[:] # ϵ-μ # -k
        end
    end
    G
end

###################################################################################################
###################################################################################################
###################################################################################################

import LatticeQM.Eigen
import LatticeQM.Utils

function green!(G::AbstractHops, H, ks::AbstractMatrix{Float64}, μ::Float64=0.0;
                T::Float64=0.01,
                multimode::Symbol=:auto,
                executor::Union{Nothing,Parallel.Executor}=nothing,
                kwargs...)
    L = size(ks, 2)
    for δL in keys(G)
        G[δL] .= 0
    end

    exec = executor === nothing ? Parallel.to_executor(multimode) : executor
    Parallel.configure_blas!(exec; verbose=false)

    energy_lock = Threads.SpinLock()
    total_energy = Ref(0.0)

    Parallel.kspace_reduce!(G, ks, exec) do local_G, _scratch, _j, k
        ϵs, U = Eigen.geteigen(H(k); kwargs...)
        green!(local_G, k, real.(ϵs) .- μ, U; T=T)
        e = Utils.groundstate_sumk(real.(ϵs), μ)
        Base.@lock energy_lock total_energy[] += e
    end

    for δL in keys(G)
        G[δL] ./= L
    end

    total_energy[] / L
end

