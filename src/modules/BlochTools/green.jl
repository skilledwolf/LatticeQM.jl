
function G!(G0::AbstractMatrix, h::Function, ks::AbstractMatrix{Float64}, mu::Float64)
    G0[:] .= zero(G0)

    for F in energies_wfs(h, ks) # @todo: this should be paralellized

        for (ϵ, ψ) in zip(F.values, eachcol(F.vectors))
            if ϵ <= μ
                G0[:] .+= ψ * ψ'
            end
        end
    end

    nothing
end

function G(h::Function, ks::AbstractMatrix{Float64}, mu::Float64)
    G0 = similar(h(ks[:,1]))
    G!(G0, h, ks, mu)

    G0
end
