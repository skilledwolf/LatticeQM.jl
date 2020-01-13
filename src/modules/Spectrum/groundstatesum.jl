@fastmath function groundstate_sumk(ϵs_k::AbstractVector{Float64}, μ::Float64=0.0)
    tmp = 0.0
    for ϵ in ϵs_k
        if ϵ <= μ
            tmp += ϵ
        end
    end

    tmp
end

function groundstate_energy(ϵs::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; kwargs...)
    # Σ = ϵs(hamiltonian; format=format)
    L = size(ks)[2]

    ϵGS = @distributed (+) for j=1:L
        groundstate_sumk(ϵs(ks[:,j]), μ)
    end

    ϵGS / L
end
