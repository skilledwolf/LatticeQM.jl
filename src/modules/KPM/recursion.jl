using Distributed
using SharedArrays
using SparseArrays
using LinearAlgebra

function kpm_step!(α2::T2, α1::T2, α0::T2, H::T1) where {T1<:AbstractMatrix{<:Complex}, T2<:AbstractVector{<:Complex}}
"""
    Iterative step in the KPM for Hermitian matrix H (with bounded spectrum [-1,1])
    to compute the next vector (and hence the next Chebyshev expansion coefficient).
"""
    α2[:] = (2.0 .* (H * α1) - α0)[:]
    α0[:] = α1[:]
    α1[:] = α2[:]
    nothing
end


function get_expansion!(μ, α::T2, H::T1) where {T1<:AbstractMatrix{<:Complex}, T2<:AbstractVector{<:Complex}}
    n = length(μ)

    α0 = deepcopy(α)
    α1 = H * α0
    αTMP = similar(α)

    μ0 = α' * α0
    μ1 = α' * α1

    μ[1] += μ0
    μ[2] += μ1

    for j_=3:2:n
        kpm_step!(αTMP, α1, α0, H) # α1 = α_{N+1},  α0 = α_{N}
        μ[j_] += 2.0 * α0' * α0 - μ0
        if j_ < n
            μ[j_+1] += 2.0 * α1' * α0 - μ1
        end
    end

    nothing
end

function get_expansion!(μ, β::T2, α::T2, H::T1) where {T1<:AbstractMatrix{<:Complex}, T2<:AbstractVector{<:Complex}}
    n = length(μ)

    α0 = deepcopy(α)
    α1 = H * α0
    αTMP = similar(α)

    μ[1] += β' * α0
    μ[2] += β' * α1

    for j_=3:n
        kpm_step!(αTMP, α1, α0, H)
        μ[j_] += β' * α1
    end

    nothing
end

function get_expansion(n::Int, args...; kwargs...)
"""
    Evaluate the KPM expansion coefficients Eq. (28) in [1] using
    the recursion for Chebyshev polynomials (see Sec. II.B).

    H is the Hamiltonian matrix, n is the expansion order, R is the number of random
    vectors to use in the stochastic trace.
    A is the operator for which we want the expectation value and will
    depenend on the physical problem (e.g. A=id for density of states).

    [1] Weisse et al, Rev. Mod. Phys. 78 275
"""
    μ = zeros(ComplexF64, n)
    get_expansion!(μ, args...; kwargs...)

    μ
end

function get_expansion_stochastic(A::T, H::T, n::Int, R::Int; mode=:normal, kernel=:jackson, kernelparams...) where T<:AbstractMatrix
"""
    Evaluate the KPM expansion coefficients Eq. (29) in [1] using
    a stochastic trace (see Sec. II.B).

    H is the Hamiltonian, n is the expansion order, R is the number of random
    vectors to use in the stochastic trace.
    A is the operator for which we want the expectation value and will
    depenend on the physical problem (e.g. A=id for density of states).

    [1] Weisse et al, Rev. Mod. Phys. 78 275
"""
#     μ = convert(SharedArray, zeros(ComplexF64, n))
    D = size(H,2)
    # αr = zeros(D)

    if mode==:normal
        myrand = randn
    elseif mode==:uniform
        myrand(args...; kwargs...) = √3.0 .* (2.0 .* rand(args...; kwargs...) .- 1.0)
    end

    μ = @distributed (+) for r_=1:R
        μ0 = zeros(ComplexF64, n)
        αR = Vector{ComplexF64}(myrand(D))
        get_expansion!(μ0, αR, H)
        μ0
    end
#     @sync @distributed for r_=1:R
#         αR = Vector{ComplexF64}(myrand(D))
#         get_expansion!(μ, A'*αR, αR, H)
#     end

    if kernel==:jackson
        jackson!(μ)
    elseif kernel==:lorentz
        lorentz!(μ)
    end

    Array(μ ./ R)
end

function get_expansion_stochastic(H::T, n::Int, R::Int; mode=:normal, kernel=:jackson, kernelparams...) where T<:AbstractMatrix
"""
    Evaluate the KPM expansion coefficients Eq. (29) in [1] using
    a stochastic trace (see Sec. II.B).

    H is the Hamiltonian, n is the expansion order, R is the number of random
    vectors to use in the stochastic trace.
    A is the operator for which we want the expectation value and will
    depenend on the physical problem (e.g. A=id for density of states).

    [1] Weisse et al, Rev. Mod. Phys. 78 275
"""
    μ = convert(SharedArray, zeros(ComplexF64, n))
    D = size(H,2)
    # αr = zeros(D)

    if mode==:normal
        myrand = randn
    elseif mode==:uniform
        myrand(args...; kwargs...) = √3.0 .* (2.0 .* rand(args...; kwargs...) .- 1.0)
    end

    μ = @distributed (+) for r_=1:R
        μ0 = zeros(ComplexF64, n)
        αR = Vector{ComplexF64}(myrand(D))
        get_expansion!(μ0, αR, H)
        μ0
    end
#     @sync @distributed for r_=1:R
#         αR = Vector{ComplexF64}(myrand(D))
#         get_expansion!(μ, αR, H)
#     end

    if kernel==:jackson
        jackson!(μ)
    elseif kernel==:lorentz
        lorentz!(μ)
    end

    Array(μ ./ R)
end
