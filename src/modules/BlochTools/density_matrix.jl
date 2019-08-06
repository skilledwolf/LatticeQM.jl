
using Distributed


function ρ_k!(ρ0::AbstractMatrix{ComplexF64}, spectrum_k, μ::Float64)
    ϵs, U = spectrum_k
    for (ϵ, ψ) in zip(ϵs, eachcol(U))
        if ϵ <= μ
            ρ0[:] .+= (ψ * ψ')[:]
        end
    end

    nothing
end

function ρ!(ρ0::AbstractMatrix{ComplexF64}, spectrum::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0)
    ρ0[:] .= zero(ρ0)[:]
    L = size(ks)[2]

    ρ0[:] = @distributed (+) for j=1:L
        ρ0 = zero(ρ0)    ## <-- it annoys me that I don't know how to get around this allocation
        ρ_k!(ρ0, spectrum(ks[:,j]), μ) #; format=format)

        ρ0
    end

    ρ0[:] .= ρ0[:] ./ L

    nothing
end

function ρ(spectrum::Function, d::Int, ks::AbstractMatrix{Float64}, μ::Float64; kwargs...)
    ρ0 = zeros(ComplexF64, d, d)
    ρ!(ρ0, spectrum, ks, μ; kwargs...)

    ρ0
end

# function ρ!(ρ0::AbstractMatrix{ComplexF64}, hamiltonian::Function, ks::AbstractMatrix{Float64}, μ::Float64=0.0; format=:auto)
#     ρ0[:] .= zero(ρ0)[:]
#
#     # if format==:auto # Decide if the matrix is dense or sparse
#     #     format = issparse(hamiltonian(ks[:,1])) ? :sparse : :dense
#     # end
#
#     L = size(ks)[2]
#
#     @inbounds for j=1:L # @todo: this should be paralellized
#         ρ_k!(ρ0, hamiltonian, ks[:,j], μ)#; format=format)
#     end
#
#     ρ0[:] .= ρ0[:] ./ L
#
#     nothing
# end
